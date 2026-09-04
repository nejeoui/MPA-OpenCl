#include "mpa_ref.h"
#include <time.h>

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static double median(double *v, int n)
{
    qsort(v, (size_t)n, sizeof(double), cmp_double);
    return (n & 1) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

static double minimum(const double *v, int n)
{
    double m = v[0];
    for (int i = 1; i < n; i++) if (v[i] < m) m = v[i];
    return m;
}

typedef struct {
    const char *name;
    const char *src;
    const char *flags;
    int         interleaved;
} Config;

#define OPT "mpaKernel_32bits_opt.cl"

static const Config CONFIGS[] = {
    { "baseline (fixed)",  "mpaKernel_32bits.cl", "",                  0 },
    { "opt (no flags)",    OPT,                   "",                  0 },
    { "only mul_hi",       OPT, "-DMPA_MULHI=1",                       0 },
    { "only reg accum",    OPT, "-DMPA_REGACC=1",                      0 },
    { "only fused CIOS",   OPT, "-DMPA_FUSED_CIOS=1",                  0 },
    { "only unroll",       OPT, "-DMPA_UNROLL=1",                      0 },
    { "only interleaved",  OPT, "-DMPA_INTERLEAVED=1",                 1 },
    { "cios+reg",          OPT, "-DMPA_FUSED_CIOS=1 -DMPA_REGACC=1",   0 },
    { "cios+reg+unroll",   OPT, "-DMPA_FUSED_CIOS=1 -DMPA_REGACC=1 -DMPA_UNROLL=1", 0 },
    { "all but mul_hi",    OPT, "-DMPA_FUSED_CIOS=1 -DMPA_REGACC=1 -DMPA_UNROLL=1 -DMPA_INTERLEAVED=1", 1 },
    { "all incl. mul_hi",  OPT, "-DMPA_MULHI=1 -DMPA_FUSED_CIOS=1 -DMPA_REGACC=1 -DMPA_UNROLL=1 -DMPA_INTERLEAVED=1", 1 },
};
#define NCONFIGS ((int)(sizeof(CONFIGS)/sizeof(CONFIGS[0])))

static size_t addr(int interleaved, size_t it, int w, int nwords, size_t items)
{
    return interleaved ? (size_t)w * items + it : it * (size_t)nwords + (size_t)w;
}

static void putOperand(const mpz_t z, uint32_t *buf, size_t it, int nwords,
                       int interleaved, size_t items)
{
    mpz_t t;
    mpz_init_set(t, z);
    for (int i = nwords - 1; i >= 0; i--) {
        buf[addr(interleaved, it, i, nwords, items)] =
            (uint32_t)(mpz_get_ui(t) & 0xFFFFFFFFu);
        mpz_tdiv_q_2exp(t, t, 32);
    }
    mpz_clear(t);
}

typedef struct {
    double kernel_s;
    double e2e_s;
    double noise;
    long   mismatches;
    int    built;
} Result;

static Result runConfig(cl_context ctx, cl_command_queue q, cl_device_id dev,
                        const Config *cfg, const Op *op, int T,
                        const mpz_t p, unsigned long m_prime,
                        size_t items, int reps, int verbose)
{
    Result R = { 0, 0, 0, 0, 0 };
    const int outWords = op->wide ? 2 * T : T;
    const size_t inBytes  = items * (size_t)T * 4;
    const size_t outBytes = items * (size_t)outWords * 4;

    uint32_t *hA   = calloc(items * (size_t)T, 4);
    uint32_t *hB   = calloc(items * (size_t)T, 4);
    uint32_t *hC   = calloc(items * (size_t)outWords, 4);
    uint32_t *hExp = calloc(items * (size_t)outWords, 4);
    uint32_t *hP   = calloc((size_t)T * 2, 4);
    cl_int opBuf[4];

    mpz_t a, b, e;
    mpz_inits(a, b, e, NULL);
    for (size_t j = 0; j < items; j++) {
        makeCase((int)(j % 4096), op->op, op->modular, p, T * 32, &MODULI[0], a, b);
        putOperand(a, hA, j, T, cfg->interleaved, items);
        putOperand(b, hB, j, T, cfg->interleaved, items);
        reference(op->op, a, b, p, T * 32, e);
        putOperand(e, hExp, j, outWords, cfg->interleaved, items);
    }
    {
        mpz_t t, r2;
        mpz_inits(t, r2, NULL);
        mpz_set(t, p);
        for (int i = T - 1; i >= 0; i--) {
            hP[i] = (uint32_t)(mpz_get_ui(t) & 0xFFFFFFFFu);
            mpz_tdiv_q_2exp(t, t, 32);
        }
        computeR2(r2, p, T * 32);
        for (int i = T - 1; i >= 0; i--) {
            hP[T + i] = (uint32_t)(mpz_get_ui(r2) & 0xFFFFFFFFu);
            mpz_tdiv_q_2exp(r2, r2, 32);
        }
        mpz_clears(t, r2, NULL);
    }
    mpz_clears(a, b, e, NULL);

    opBuf[0] = op->op; opBuf[1] = 32; opBuf[2] = T * 32;
    opBuf[3] = (cl_int)(uint32_t)m_prime;

    size_t srcLen;
    char *src = readFile(cfg->src, &srcLen);
    char opts[512];
    snprintf(opts, sizeof(opts), "-I. -DWORDLENGTH_T=%d %s", T, cfg->flags);

    cl_int err;
    cl_program prog = clCreateProgramWithSource(ctx, 1, (const char **)&src, &srcLen, &err);
    CHECK(err);
    if (clBuildProgram(prog, 1, &dev, opts, NULL, NULL) != CL_SUCCESS) {
        size_t ln = 0;
        clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &ln);
        char *log = malloc(ln + 1);
        clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, ln, log, NULL);
        log[ln] = 0;
        fprintf(stderr, "build failed [%s] %s:\n%s\n", cfg->name, opts, log);
        free(log); free(src);
        goto cleanup_host;
    }
    free(src);
    cl_kernel kern = clCreateKernel(prog, "mpaKernel", &err); CHECK(err);
    R.built = 1;

    cl_mem dA = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  inBytes,  NULL, &err); CHECK(err);
    cl_mem dB = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  inBytes,  NULL, &err); CHECK(err);
    cl_mem dC = clCreateBuffer(ctx, CL_MEM_READ_WRITE, outBytes, NULL, &err); CHECK(err);
    cl_mem dO = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  sizeof(opBuf), NULL, &err); CHECK(err);
    cl_mem dP = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  (size_t)T*2*4, NULL, &err); CHECK(err);

    CHECK(clSetKernelArg(kern, 0, sizeof(cl_mem), &dA));
    CHECK(clSetKernelArg(kern, 1, sizeof(cl_mem), &dB));
    CHECK(clSetKernelArg(kern, 2, sizeof(cl_mem), &dC));
    CHECK(clSetKernelArg(kern, 3, sizeof(cl_mem), &dO));
    CHECK(clSetKernelArg(kern, 4, sizeof(cl_mem), &dP));

    CHECK(clEnqueueWriteBuffer(q, dO, CL_TRUE, 0, sizeof(opBuf), opBuf, 0, NULL, NULL));
    CHECK(clEnqueueWriteBuffer(q, dP, CL_TRUE, 0, (size_t)T*2*4, hP, 0, NULL, NULL));

    CHECK(clEnqueueWriteBuffer(q, dA, CL_TRUE, 0, inBytes, hA, 0, NULL, NULL));
    CHECK(clEnqueueWriteBuffer(q, dB, CL_TRUE, 0, inBytes, hB, 0, NULL, NULL));
    CHECK(clEnqueueNDRangeKernel(q, kern, 1, NULL, &items, NULL, 0, NULL, NULL));
    CHECK(clFinish(q));
    CHECK(clEnqueueReadBuffer(q, dC, CL_TRUE, 0, outBytes, hC, 0, NULL, NULL));

    for (size_t j = 0; j < items; j++)
        for (int w = 0; w < outWords; w++) {
            size_t k = addr(cfg->interleaved, j, w, outWords, items);
            if (hC[k] != hExp[k]) { R.mismatches++; break; }
        }

    if (R.mismatches) {
        if (verbose)
            fprintf(stderr, "  [%s] %ld/%zu wrong - not timed\n",
                    cfg->name, R.mismatches, items);
        goto cleanup_dev;
    }

    for (int r = 0; r < 2; r++) {
        CHECK(clEnqueueNDRangeKernel(q, kern, 1, NULL, &items, NULL, 0, NULL, NULL));
        CHECK(clFinish(q));
    }

    double *tk = malloc(sizeof(double) * (size_t)reps);
    double *te = malloc(sizeof(double) * (size_t)reps);
    for (int r = 0; r < reps; r++) {
        double t0 = now_s();
        CHECK(clEnqueueNDRangeKernel(q, kern, 1, NULL, &items, NULL, 0, NULL, NULL));
        CHECK(clFinish(q));
        tk[r] = now_s() - t0;

        double t1 = now_s();
        CHECK(clEnqueueWriteBuffer(q, dA, CL_FALSE, 0, inBytes, hA, 0, NULL, NULL));
        CHECK(clEnqueueWriteBuffer(q, dB, CL_FALSE, 0, inBytes, hB, 0, NULL, NULL));
        CHECK(clEnqueueNDRangeKernel(q, kern, 1, NULL, &items, NULL, 0, NULL, NULL));
        CHECK(clEnqueueReadBuffer(q, dC, CL_TRUE, 0, outBytes, hC, 0, NULL, NULL));
        CHECK(clFinish(q));
        te[r] = now_s() - t1;
    }
    R.kernel_s = minimum(tk, reps);
    R.e2e_s    = minimum(te, reps);
    {
        double *cp = malloc(sizeof(double) * (size_t)reps);
        memcpy(cp, tk, sizeof(double) * (size_t)reps);
        double med = median(cp, reps);
        R.noise = (R.kernel_s > 0) ? (med - R.kernel_s) / R.kernel_s : 0.0;
        free(cp);
    }
    free(tk); free(te);

cleanup_dev:
    clReleaseMemObject(dA); clReleaseMemObject(dB); clReleaseMemObject(dC);
    clReleaseMemObject(dO); clReleaseMemObject(dP);
    clReleaseKernel(kern);
    clReleaseProgram(prog);
cleanup_host:
    free(hA); free(hB); free(hC); free(hExp); free(hP);
    return R;
}

int main(int argc, char **argv)
{
    size_t items = 16384;
    int reps = 15, bits = 256, allOps = 0;
    const char *onlyOp = NULL, *csvPath = NULL;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--items") && i+1 < argc) items = strtoul(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--reps")  && i+1 < argc) reps  = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--bits")  && i+1 < argc) bits  = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--op")    && i+1 < argc) onlyOp = argv[++i];
        else if (!strcmp(argv[i], "--csv")   && i+1 < argc) csvPath = argv[++i];
        else if (!strcmp(argv[i], "--all-ops")) allOps = 1;
        else {
            fprintf(stderr, "usage: %s [--items N] [--reps R] [--bits 256|512|1024|2048]\n"
                            "          [--op NAME] [--all-ops] [--csv FILE]\n", argv[0]);
            return 2;
        }
    }
    if (reps < 1) reps = 1;
    rngState = 20260902;

    const Modulus *mod = NULL;
    for (int m = 0; m < NMODULI; m++)
        if (MODULI[m].bits == bits && !MODULI[m].facA) { mod = &MODULI[m]; break; }
    if (!mod) { fprintf(stderr, "no prime modulus of %d bits in the table\n", bits); return 2; }
    const int T = bits / 32;

    mpz_t p, base, inv, mp;
    mpz_inits(p, base, inv, mp, NULL);
    mpz_set_str(p, mod->hex, 16);
    mpz_ui_pow_ui(base, 2, 32);
    mpz_invert(inv, p, base);
    mpz_sub(mp, base, inv);
    unsigned long m_prime = mpz_get_ui(mp);

    cl_platform_id plat; cl_device_id dev;
    pickDevice(&plat, &dev);
    char nameBuf[256] = {0}, verBuf[128] = {0};
    cl_uint cus = 0;
    cl_device_type dtype = 0;
    clGetDeviceInfo(dev, CL_DEVICE_NAME, sizeof(nameBuf), nameBuf, NULL);
    clGetDeviceInfo(dev, CL_DEVICE_VERSION, sizeof(verBuf), verBuf, NULL);
    clGetDeviceInfo(dev, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cus), &cus, NULL);
    clGetDeviceInfo(dev, CL_DEVICE_TYPE, sizeof(dtype), &dtype, NULL);

    cl_int err;
    cl_context ctx = clCreateContext(NULL, 1, &dev, NULL, NULL, &err); CHECK(err);
    cl_command_queue q = clCreateCommandQueue(ctx, dev, 0, &err);      CHECK(err);

    printf("device  : [%s] %s (%u compute units)\n", devTypeName(dtype), nameBuf, cus);
    printf("version : %s\n", verBuf);
    printf("modulus : %s (%d bits, T=%d words of 32 bits)\n", mod->name, bits, T);
    printf("workload: %zu items, %d timed repetitions, MINIMUM reported\n", items, reps);
    printf("          'noise' is (median-min)/min: high values mean a busy host\n\n");

    FILE *csv = NULL;
    if (csvPath) {
        csv = fopen(csvPath, "w");
        if (csv) fprintf(csv, "operation,config,items,bits,kernel_s,e2e_s,noise,kernel_ops_s,e2e_ops_s,speedup\n");
    }

    const char *defaultOps[] = { "ADDMOD", "MULTIPLYOPERANDSCANNING",
                                 "MULTIPLYPRODUCTSCANNING", "MONTGOMERYMULTIPLICATION" };
    int broken = 0;

    for (int o = 0; o < NOPS; o++) {
        const Op *op = &OPS[o];
        if (onlyOp) {
            if (strcmp(op->name, onlyOp)) continue;
        } else if (!allOps) {
            int want = 0;
            for (int d = 0; d < 4; d++) if (!strcmp(op->name, defaultOps[d])) want = 1;
            if (!want) continue;
        }

        printf("%s\n", op->name);
        printf("  %-18s %11s %8s %11s %13s %13s %9s\n",
               "config", "kernel ms", "noise", "e2e ms", "kernel op/s", "e2e op/s", "vs base");
        printf("  %s\n", "-------------------------------------------------------------------------------------");

        double baseKernel = 0;
        for (int c = 0; c < NCONFIGS; c++) {
            if (op->ext && c == 0) continue;
            Result R = runConfig(ctx, q, dev, &CONFIGS[c], op, T, p, m_prime, items, reps, 1);
            if (!R.built) { printf("  %-18s %s\n", CONFIGS[c].name, "BUILD FAILED"); broken++; continue; }
            if (R.mismatches) {
                printf("  %-18s %s%ld/%zu results WRONG - not timed%s\n",
                       CONFIGS[c].name, RED, R.mismatches, items, OFF);
                broken++;
                continue;
            }
            double kops = (double)items / R.kernel_s;
            double eops = (double)items / R.e2e_s;
            if (baseKernel == 0) baseKernel = R.kernel_s;
            printf("  %-18s %11.3f %7.0f%% %11.3f %13.3e %13.3e %8.2fx\n",
                   CONFIGS[c].name, R.kernel_s * 1e3, R.noise * 100.0, R.e2e_s * 1e3,
                   kops, eops, baseKernel > 0 ? baseKernel / R.kernel_s : 1.0);
            if (csv)
                fprintf(csv, "%s,\"%s\",%zu,%d,%.9f,%.9f,%.4f,%.3f,%.3f,%.4f\n",
                        op->name, CONFIGS[c].name, items, bits, R.kernel_s, R.e2e_s,
                        R.noise, kops, eops, baseKernel > 0 ? baseKernel / R.kernel_s : 1.0);
        }
        printf("\n");
    }

    if (csv) fclose(csv);
    clReleaseCommandQueue(q);
    clReleaseContext(ctx);
    mpz_clears(p, base, inv, mp, NULL);

    if (broken) {
        printf("%s%d configuration(s) failed to build or verify%s\n", RED, broken, OFF);
        return 1;
    }
    return 0;
}
