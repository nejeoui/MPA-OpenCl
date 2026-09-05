#include "mpa_ref.h"
#include <openssl/bn.h>
#include <openssl/crypto.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <ctype.h>
#include <signal.h>
#include <time.h>
#ifdef __APPLE__
#include <sys/sysctl.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    char   gpuName[256], gpuVendor[160], clVersion[160], driver[160];
    cl_ulong gmem, lmem, cache, maxAlloc;
    cl_uint  cus, clockMHz;
    size_t   maxWG;
    char   cpu[256], os[256], kernel[128], arch[64];
    long   cores;
    double ramGB;
    int    threads;
} Env;

typedef struct {
    long   items;
    double gmp1, gmpN, ossl;
} CpuRow;

typedef struct {
    int    attempted, built;
    long   items, mismatches;
    double kernel_s, e2e_s;
} GpuCell;

typedef struct {
    char   mod[64], op[48];
    double seconds;
    long   items;
} CgbnRow;

static CpuRow  g_cpu[NMODULI][NOPS];
static GpuCell g_gpu[NVARIANTS][NMODULI][NOPS];
static CgbnRow g_cgbn[256];
static int     g_ncgbn = 0;
static Env     g_env;
static char    g_reportPath[512], g_csvPath[512];
static volatile sig_atomic_t g_interrupted = 0;

static void onSigint(int s) { (void)s; g_interrupted = 1; }

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static double minimum(const double *v, int n)
{
    double m = v[0];
    for (int i = 1; i < n; i++) if (v[i] < m) m = v[i];
    return m;
}

static void sanitize(const char *in, char *out, size_t n)
{
    size_t j = 0;
    int us = 0;
    for (size_t i = 0; in[i] && j + 1 < n; i++) {
        unsigned char c = (unsigned char)in[i];
        if (isalnum(c) || c == '.' || c == '-') { out[j++] = (char)c; us = 0; }
        else if (!us && j) { out[j++] = '_'; us = 1; }
    }
    while (j && out[j-1] == '_') j--;
    out[j] = 0;
    if (!j) snprintf(out, n, "UnknownGPU");
}

static void detectCpu(Env *e)
{
    snprintf(e->cpu, sizeof e->cpu, "unknown");
    e->ramGB = 0;
#ifdef __APPLE__
    size_t sz = sizeof e->cpu;
    sysctlbyname("machdep.cpu.brand_string", e->cpu, &sz, NULL, 0);
    uint64_t mem = 0; sz = sizeof mem;
    if (sysctlbyname("hw.memsize", &mem, &sz, NULL, 0) == 0)
        e->ramGB = (double)mem / 1073741824.0;
#else
    FILE *f = fopen("/proc/cpuinfo", "r");
    if (f) {
        char line[512];
        while (fgets(line, sizeof line, f)) {
            char *c = strchr(line, ':');
            if (!c) continue;
            if (!strncmp(line, "model name", 10) || !strncmp(line, "Model name", 10) ||
                !strncmp(line, "Hardware", 8)) {
                c++; while (*c == ' ' || *c == '\t') c++;
                c[strcspn(c, "\n")] = 0;
                snprintf(e->cpu, sizeof e->cpu, "%s", c);
                break;
            }
        }
        fclose(f);
    }
    f = fopen("/proc/meminfo", "r");
    if (f) {
        char line[256];
        while (fgets(line, sizeof line, f)) {
            unsigned long kb;
            if (sscanf(line, "MemTotal: %lu kB", &kb) == 1) {
                e->ramGB = (double)kb / 1048576.0;
                break;
            }
        }
        fclose(f);
    }
#endif
    e->cores = sysconf(_SC_NPROCESSORS_ONLN);
}

static void detectOs(Env *e)
{
    struct utsname u;
    e->os[0] = e->kernel[0] = e->arch[0] = 0;
    FILE *f = fopen("/etc/os-release", "r");
    if (f) {
        char line[512];
        while (fgets(line, sizeof line, f)) {
            if (!strncmp(line, "PRETTY_NAME=", 12)) {
                char *v = line + 12;
                if (*v == '"') v++;
                char *q = strrchr(v, '"');
                if (q) *q = 0; else v[strcspn(v, "\n")] = 0;
                snprintf(e->os, sizeof e->os, "%s", v);
                break;
            }
        }
        fclose(f);
    }
    if (uname(&u) == 0) {
        snprintf(e->kernel, sizeof e->kernel, "%s", u.release);
        snprintf(e->arch, sizeof e->arch, "%s", u.machine);
        if (!e->os[0]) snprintf(e->os, sizeof e->os, "%s %s", u.sysname, u.release);
    }
    if (!e->os[0]) snprintf(e->os, sizeof e->os, "unknown");
}

static void detectGpu(Env *e, cl_device_id d)
{
    clGetDeviceInfo(d, CL_DEVICE_NAME,                sizeof e->gpuName,   e->gpuName,   NULL);
    clGetDeviceInfo(d, CL_DEVICE_VENDOR,              sizeof e->gpuVendor, e->gpuVendor, NULL);
    clGetDeviceInfo(d, CL_DEVICE_VERSION,             sizeof e->clVersion, e->clVersion, NULL);
    clGetDeviceInfo(d, CL_DRIVER_VERSION,             sizeof e->driver,    e->driver,    NULL);
    clGetDeviceInfo(d, CL_DEVICE_GLOBAL_MEM_SIZE,      sizeof e->gmem,     &e->gmem,     NULL);
    clGetDeviceInfo(d, CL_DEVICE_LOCAL_MEM_SIZE,       sizeof e->lmem,     &e->lmem,     NULL);
    clGetDeviceInfo(d, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,sizeof e->cache,    &e->cache,    NULL);
    clGetDeviceInfo(d, CL_DEVICE_MAX_MEM_ALLOC_SIZE,   sizeof e->maxAlloc, &e->maxAlloc, NULL);
    clGetDeviceInfo(d, CL_DEVICE_MAX_COMPUTE_UNITS,    sizeof e->cus,      &e->cus,      NULL);
    clGetDeviceInfo(d, CL_DEVICE_MAX_CLOCK_FREQUENCY,  sizeof e->clockMHz, &e->clockMHz, NULL);
    clGetDeviceInfo(d, CL_DEVICE_MAX_WORK_GROUP_SIZE,  sizeof e->maxWG,    &e->maxWG,    NULL);
}

typedef struct { mpz_t lim, Rinv, tmp; int bits; } GmpCtx;

static void gmpOp(int op, mpz_t r, const mpz_t a, const mpz_t b,
                  const mpz_t p, GmpCtx *c)
{
    switch (op) {
    case ADD:            mpz_add(r, a, b); mpz_fdiv_r_2exp(r, r, (mp_bitcnt_t)c->bits); break;
    case SUBTRACT:       mpz_sub(r, a, b); mpz_fdiv_r_2exp(r, r, (mp_bitcnt_t)c->bits); break;
    case ADDMOD:         mpz_add(r, a, b); mpz_mod(r, r, p); break;
    case SUBTRACTMOD:    mpz_sub(r, a, b); mpz_mod(r, r, p); break;
    case MULTIPLYOPERANDSCANNING:
    case MULTIPLYPRODUCTSCANNING: mpz_mul(r, a, b); break;
    case MONTGOMERYMULTIPLICATION:
        mpz_mul(r, a, b); mpz_mul(r, r, c->Rinv); mpz_mod(r, r, p); break;
    case COMPARE:        mpz_set_si(r, mpz_cmp(a, b)); break;
    case REDUCE:         mpz_mod(r, a, p); break;
    case MODMUL:
    case MODMUL_R2:      mpz_mul(r, a, b); mpz_mod(r, r, p); break;
    case MODEXP:         mpz_powm(r, a, b, p); break;
    case EXPONENTIATION: mpz_powm(r, a, b, c->lim); break;
    case DIVIDE:
        if (mpz_sgn(b)) mpz_tdiv_qr(r, c->tmp, a, b); else mpz_set_ui(r, 0);
        break;
    case ISQRT:          mpz_sqrt(r, a); break;
    default:             mpz_set_ui(r, 0); break;
    }
}

static int osslSupports(int op)
{
    switch (op) {
    case ISQRT: return 0;
    default:    return 1;
    }
}

static void osslOp(int op, BIGNUM *r, BIGNUM *scratch, const BIGNUM *a,
                   const BIGNUM *b, const BIGNUM *p, const BIGNUM *lim,
                   BN_MONT_CTX *mont, BN_CTX *ctx, int bits)
{
    switch (op) {
    case ADD:            BN_add(r, a, b); BN_mask_bits(r, bits); break;
    case SUBTRACT:       BN_sub(r, a, b); BN_mask_bits(r, bits); break;
    case ADDMOD:         BN_mod_add(r, a, b, p, ctx); break;
    case SUBTRACTMOD:    BN_mod_sub(r, a, b, p, ctx); break;
    case MULTIPLYOPERANDSCANNING:
    case MULTIPLYPRODUCTSCANNING: BN_mul(r, a, b, ctx); break;
    case MONTGOMERYMULTIPLICATION:
        if (mont) BN_mod_mul_montgomery(r, a, b, mont, ctx); break;
    case COMPARE:        BN_set_word(r, (BN_ULONG)(BN_cmp(a, b) + 1)); break;
    case REDUCE:         BN_nnmod(r, a, p, ctx); break;
    case MODMUL:
    case MODMUL_R2:      BN_mod_mul(r, a, b, p, ctx); break;
    case MODEXP:         BN_mod_exp(r, a, b, p, ctx); break;
    case EXPONENTIATION: BN_mod_exp(r, a, b, lim, ctx); break;
    case DIVIDE:
        if (!BN_is_zero(b)) BN_div(r, scratch, a, b, ctx); break;
    default: break;
    }
}

static const char *modShort(const Modulus *m) { return m->name; }

static void writeReport(int items, int reps, double elapsed, int complete);

int main(int argc, char **argv)
{
    int items = 20000, reps = 5, budget = 0;
    const char *only = NULL;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--items")  && i+1 < argc) items  = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--reps")   && i+1 < argc) reps   = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--budget") && i+1 < argc) budget = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--variant")&& i+1 < argc) only   = argv[++i];
        else if (!strcmp(argv[i], "--verbose")) g_verbose = 1;
        else {
            fprintf(stderr,
                "usage: %s [--items N] [--reps N] [--budget SECONDS] "
                "[--variant w8|w16|w32|w32-opt|w32-o64] [--verbose]\n", argv[0]);
            return 2;
        }
    }
    if (items < 64) items = 64;
    if (reps  < 1)  reps  = 1;

    signal(SIGINT, onSigint);

    rngState = 88172645463325252ULL;

    cl_platform_id plat; cl_device_id dev;
    pickDevice(&plat, &dev);

    detectGpu(&g_env, dev);
    detectCpu(&g_env);
    detectOs(&g_env);
#ifdef _OPENMP
    g_env.threads = omp_get_max_threads();
#else
    g_env.threads = 1;
#endif

    char safe[256];
    sanitize(g_env.gpuName, safe, sizeof safe);
    snprintf(g_reportPath, sizeof g_reportPath, "%s_Report.md", safe);
    snprintf(g_csvPath,    sizeof g_csvPath,    "%s_Report.csv", safe);

    fprintf(stderr, "GPU     : %s (%s)\n", g_env.gpuName, g_env.gpuVendor);
    fprintf(stderr, "CPU     : %s, %ld cores, %.1f GB\n", g_env.cpu, g_env.cores, g_env.ramGB);
    fprintf(stderr, "OS      : %s\n", g_env.os);
    fprintf(stderr, "report  : %s\n\n", g_reportPath);

    {
        FILE *f = fopen("cgbn_results.tsv", "r");
        if (f) {
            char line[512];
            while (g_ncgbn < 256 && fgets(line, sizeof line, f)) {
                if (line[0] == '#' || line[0] == '\n') continue;
                CgbnRow r;
                if (sscanf(line, "%63s %47s %ld %lf", r.mod, r.op, &r.items, &r.seconds) == 4)
                    g_cgbn[g_ncgbn++] = r;
            }
            fclose(f);
            fprintf(stderr, "cgbn    : %d rows from cgbn_results.tsv\n\n", g_ncgbn);
        } else {
            fprintf(stderr, "cgbn    : cgbn_results.tsv absent, CGBN columns will read n/a\n\n");
        }
    }

    const double t_start = now_s();

    cl_int err;
    cl_context ctx = clCreateContext(NULL, 1, &dev, NULL, NULL, &err); CHECK(err);
    cl_command_queue q = clCreateCommandQueue(ctx, dev, 0, &err); CHECK(err);

    for (int m = 0; m < NMODULI && !g_interrupted; m++) {
        const Modulus *mod = &MODULI[m];
        mpz_t p; mpz_init(p); mpz_set_str(p, mod->hex, 16);

        GmpCtx gc; mpz_inits(gc.lim, gc.Rinv, gc.tmp, NULL);
        gc.bits = mod->bits;
        mpz_ui_pow_ui(gc.lim, 2, (unsigned long)mod->bits);
        if (mpz_invert(gc.Rinv, gc.lim, p) == 0) mpz_set_ui(gc.Rinv, 1);

        BIGNUM *bp = NULL, *blim = NULL;
        BN_hex2bn(&bp, mod->hex);
        {
            char *h = mpz_get_str(NULL, 16, gc.lim);
            BN_hex2bn(&blim, h); free(h);
        }
        BN_CTX *mctx = BN_CTX_new();
        BN_MONT_CTX *mont = BN_MONT_CTX_new();
        if (!BN_MONT_CTX_set(mont, bp, mctx)) { BN_MONT_CTX_free(mont); mont = NULL; }

        for (int o = 0; o < NOPS && !g_interrupted; o++) {
            const Op *op = &OPS[o];
            int div = op->cost * (mod->bits / 256 > 1 ? mod->bits / 256 : 1);
            long n = items / (div > 1 ? div : 1);
            if (n < 64) n = 64;

            rngState = 88172645463325252ULL + (uint64_t)m * 1000u + (uint64_t)o;
            mpz_t *A = malloc(sizeof(mpz_t) * (size_t)n);
            mpz_t *B = malloc(sizeof(mpz_t) * (size_t)n);
            mpz_t *O = malloc(sizeof(mpz_t) * (size_t)n);
            for (long j = 0; j < n; j++) {
                mpz_inits(A[j], B[j], O[j], NULL);
                makeCase((int)j, op->op, op->modular, p, mod->bits, mod, A[j], B[j]);
            }

            double *t = malloc(sizeof(double) * (size_t)reps);
            for (int r = 0; r < reps; r++) {
                double t0 = now_s();
                for (long j = 0; j < n; j++) gmpOp(op->op, O[j], A[j], B[j], p, &gc);
                t[r] = now_s() - t0;
            }
            g_cpu[m][o].gmp1 = minimum(t, reps);
            g_cpu[m][o].items = n;

#ifdef _OPENMP
            for (int r = 0; r < reps; r++) {
                double t0 = now_s();
#pragma omp parallel num_threads(g_env.threads)
                {
                    GmpCtx lc; mpz_inits(lc.lim, lc.Rinv, lc.tmp, NULL);
                    lc.bits = gc.bits;
                    mpz_set(lc.lim, gc.lim); mpz_set(lc.Rinv, gc.Rinv);
                    mpz_t lr; mpz_init(lr);
#pragma omp for schedule(static)
                    for (long j = 0; j < n; j++) gmpOp(op->op, lr, A[j], B[j], p, &lc);
                    mpz_clear(lr);
                    mpz_clears(lc.lim, lc.Rinv, lc.tmp, NULL);
                }
                t[r] = now_s() - t0;
            }
            g_cpu[m][o].gmpN = minimum(t, reps);
#else
            g_cpu[m][o].gmpN = -1;
#endif

            if (osslSupports(op->op)) {
                BIGNUM **bA = malloc(sizeof(BIGNUM*) * (size_t)n);
                BIGNUM **bB = malloc(sizeof(BIGNUM*) * (size_t)n);
                for (long j = 0; j < n; j++) {
                    char *ha = mpz_get_str(NULL, 16, A[j]);
                    char *hb = mpz_get_str(NULL, 16, B[j]);
                    bA[j] = NULL; bB[j] = NULL;
                    BN_hex2bn(&bA[j], ha); BN_hex2bn(&bB[j], hb);
                    free(ha); free(hb);
                }
                for (int r = 0; r < reps; r++) {
                    double t0 = now_s();
#ifdef _OPENMP
#pragma omp parallel num_threads(g_env.threads)
#endif
                    {
                        BN_CTX *c = BN_CTX_new();
                        BIGNUM *r1 = BN_new(), *r2 = BN_new();
                        BN_MONT_CTX *lm = BN_MONT_CTX_new();
                        if (!BN_MONT_CTX_set(lm, bp, c)) { BN_MONT_CTX_free(lm); lm = NULL; }
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
                        for (long j = 0; j < n; j++)
                            osslOp(op->op, r1, r2, bA[j], bB[j], bp, blim, lm, c, mod->bits);
                        if (lm) BN_MONT_CTX_free(lm);
                        BN_free(r1); BN_free(r2); BN_CTX_free(c);
                    }
                    t[r] = now_s() - t0;
                }
                g_cpu[m][o].ossl = minimum(t, reps);
                for (long j = 0; j < n; j++) { BN_free(bA[j]); BN_free(bB[j]); }
                free(bA); free(bB);
            } else {
                g_cpu[m][o].ossl = -1;
            }

            free(t);
            for (long j = 0; j < n; j++) mpz_clears(A[j], B[j], O[j], NULL);
            free(A); free(B); free(O);

            fprintf(stderr, "  cpu  %-18s %-26s n=%-7ld gmp1=%.4fs gmpN=%.4fs ossl=%.4fs\n",
                    modShort(mod), op->name, n, g_cpu[m][o].gmp1,
                    g_cpu[m][o].gmpN, g_cpu[m][o].ossl);
        }

        if (mont) BN_MONT_CTX_free(mont);
        BN_CTX_free(mctx); BN_free(bp); BN_free(blim);
        mpz_clears(gc.lim, gc.Rinv, gc.tmp, NULL);
        mpz_clear(p);
    }

    for (int v = 0; v < NVARIANTS && !g_interrupted; v++) {
        const Variant *var = &VARIANTS[v];
        if (only && strcmp(only, var->name)) continue;

        size_t srcLen;
        char *src = readFile(var->cl, &srcLen);

        for (int m = 0; m < NMODULI && !g_interrupted; m++) {
            const Modulus *mod = &MODULI[m];
            const int T = mod->bits / var->wbits;
            const size_t wsz = (size_t)(var->wbits / 8);

            if (budget > 0 && now_s() - t_start > budget) {
                fprintf(stderr, "budget exhausted, stopping\n");
                g_interrupted = 1;
                break;
            }

            char opts[512];
            snprintf(opts, sizeof opts, "-I. -DWORDLENGTH_T=%d %s", T, var->flags);

            cl_program prog = clCreateProgramWithSource(ctx, 1, (const char**)&src, &srcLen, &err);
            if (err != CL_SUCCESS) continue;
            if (clBuildProgram(prog, 1, &dev, opts, NULL, NULL) != CL_SUCCESS) {
                size_t ln = 0;
                clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &ln);
                char *log = malloc(ln + 1);
                clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, ln, log, NULL);
                log[ln] = 0;
                fprintf(stderr, "  BUILD FAILED %s/%s: %s\n", var->name, mod->name, log);
                free(log);
                clReleaseProgram(prog);
                continue;
            }
            cl_kernel kern = clCreateKernel(prog, "mpaKernel", &err);
            if (err != CL_SUCCESS) { clReleaseProgram(prog); continue; }

            mpz_t p; mpz_init(p); mpz_set_str(p, mod->hex, 16);
            unsigned long mprime;
            {
                mpz_t base, inv, mp;
                mpz_inits(base, inv, mp, NULL);
                mpz_ui_pow_ui(base, 2, (unsigned long)var->wbits);
                if (mpz_invert(inv, p, base) == 0) mpz_set_ui(inv, 1);
                mpz_sub(mp, base, inv);
                mprime = mpz_get_ui(mp);
                mpz_clears(base, inv, mp, NULL);
            }

            for (int o = 0; o < NOPS && !g_interrupted; o++) {
                const Op *op = &OPS[o];
                GpuCell *cell = &g_gpu[v][m][o];
                if (op->ext && !var->ext) continue;
                cell->attempted = 1;

                const long n = g_cpu[m][o].items;
                if (n <= 0) { cell->attempted = 0; continue; }
                const int outWords = op->wide ? 2*T : T;
                cell->items = n;

                void *hA = calloc((size_t)n * T, wsz);
                void *hB = calloc((size_t)n * T, wsz);
                void *hC = calloc((size_t)n * outWords, wsz);
                void *hE = calloc((size_t)n * outWords, wsz);
                void *hP = calloc((size_t)T * 2, wsz);

                rngState = 88172645463325252ULL + (uint64_t)m * 1000u + (uint64_t)o;
                mpz_t a, b, e; mpz_inits(a, b, e, NULL);
                for (long j = 0; j < n; j++) {
                    makeCase((int)j, op->op, op->modular, p, mod->bits, mod, a, b);
                    mpzToWords(a, hA, (size_t)j*T, T, var->wbits);
                    mpzToWords(b, hB, (size_t)j*T, T, var->wbits);
                    reference(op->op, a, b, p, mod->bits, e);
                    mpzToWords(e, hE, (size_t)j*outWords, outWords, var->wbits);
                }
                mpzToWords(p, hP, 0, T, var->wbits);
                { mpz_t r2; mpz_init(r2); computeR2(r2, p, mod->bits);
                  mpzToWords(r2, hP, (size_t)T, T, var->wbits); mpz_clear(r2); }
                mpz_clears(a, b, e, NULL);

                cl_int ob[4] = { op->op, var->wbits, mod->bits, (cl_int)(uint32_t)mprime };
                const size_t inB = (size_t)n*T*wsz, outB = (size_t)n*outWords*wsz, pB = (size_t)T*2*wsz;

                cl_mem dA = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  inB,  NULL, &err);
                cl_mem dB = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  inB,  NULL, &err);
                cl_mem dC = clCreateBuffer(ctx, CL_MEM_READ_WRITE, outB, NULL, &err);
                cl_mem dO = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  sizeof ob, NULL, &err);
                cl_mem dP = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  pB,   NULL, &err);
                if (!dA || !dB || !dC || !dO || !dP) {
                    fprintf(stderr, "  alloc failed %s/%s/%s (%.1f MiB needed)\n",
                            var->name, mod->name, op->name,
                            (double)(2*inB + outB + pB) / 1048576.0);
                    goto cleanup;
                }

                clSetKernelArg(kern, 0, sizeof(cl_mem), &dA);
                clSetKernelArg(kern, 1, sizeof(cl_mem), &dB);
                clSetKernelArg(kern, 2, sizeof(cl_mem), &dC);
                clSetKernelArg(kern, 3, sizeof(cl_mem), &dO);
                clSetKernelArg(kern, 4, sizeof(cl_mem), &dP);

                clEnqueueWriteBuffer(q, dO, CL_TRUE, 0, sizeof ob, ob, 0, NULL, NULL);
                clEnqueueWriteBuffer(q, dP, CL_TRUE, 0, pB,  hP, 0, NULL, NULL);
                clEnqueueWriteBuffer(q, dA, CL_TRUE, 0, inB, hA, 0, NULL, NULL);
                clEnqueueWriteBuffer(q, dB, CL_TRUE, 0, inB, hB, 0, NULL, NULL);

                size_t global = (size_t)n;
                if (clEnqueueNDRangeKernel(q, kern, 1, NULL, &global, NULL, 0, NULL, NULL) != CL_SUCCESS
                    || clFinish(q) != CL_SUCCESS) {
                    fprintf(stderr, "  launch failed %s/%s/%s\n", var->name, mod->name, op->name);
                    goto cleanup;
                }
                cell->built = 1;
                clEnqueueReadBuffer(q, dC, CL_TRUE, 0, outB, hC, 0, NULL, NULL);

                for (long j = 0; j < n; j++)
                    for (int w = 0; w < outWords; w++)
                        if (loadWord(hC, (size_t)j*outWords + w, var->wbits) !=
                            loadWord(hE, (size_t)j*outWords + w, var->wbits)) { cell->mismatches++; break; }

                for (int r = 0; r < 2; r++) {
                    clEnqueueNDRangeKernel(q, kern, 1, NULL, &global, NULL, 0, NULL, NULL);
                    clFinish(q);
                }
                double *tk = malloc(sizeof(double)*(size_t)reps);
                double *te = malloc(sizeof(double)*(size_t)reps);
                for (int r = 0; r < reps; r++) {
                    double t0 = now_s();
                    clEnqueueNDRangeKernel(q, kern, 1, NULL, &global, NULL, 0, NULL, NULL);
                    clFinish(q);
                    tk[r] = now_s() - t0;

                    double t1 = now_s();
                    clEnqueueWriteBuffer(q, dA, CL_FALSE, 0, inB, hA, 0, NULL, NULL);
                    clEnqueueWriteBuffer(q, dB, CL_FALSE, 0, inB, hB, 0, NULL, NULL);
                    clEnqueueNDRangeKernel(q, kern, 1, NULL, &global, NULL, 0, NULL, NULL);
                    clEnqueueReadBuffer(q, dC, CL_TRUE, 0, outB, hC, 0, NULL, NULL);
                    clFinish(q);
                    te[r] = now_s() - t1;
                }
                cell->kernel_s = minimum(tk, reps);
                cell->e2e_s    = minimum(te, reps);
                free(tk); free(te);

                fprintf(stderr, "  gpu  %-8s %-18s %-26s n=%-7ld %.6fs %s\n",
                        var->name, modShort(mod), op->name, n, cell->kernel_s,
                        cell->mismatches ? "MISMATCH" : "ok");

cleanup:
                if (dA) clReleaseMemObject(dA); if (dB) clReleaseMemObject(dB);
                if (dC) clReleaseMemObject(dC); if (dO) clReleaseMemObject(dO);
                if (dP) clReleaseMemObject(dP);
                free(hA); free(hB); free(hC); free(hE); free(hP);
            }

            mpz_clear(p);
            clReleaseKernel(kern);
            clReleaseProgram(prog);
        }
        free(src);
        writeReport(items, reps, now_s() - t_start, 0);
    }

    writeReport(items, reps, now_s() - t_start, !g_interrupted);

    clReleaseCommandQueue(q);
    clReleaseContext(ctx);

    fprintf(stderr, "\nwrote %s and %s\n", g_reportPath, g_csvPath);
    return 0;
}

static double cgbnLookup(const char *mod, const char *op, long *items)
{
    for (int i = 0; i < g_ncgbn; i++)
        if (!strcmp(g_cgbn[i].mod, mod) && !strcmp(g_cgbn[i].op, op)) {
            if (items) *items = g_cgbn[i].items;
            return g_cgbn[i].seconds;
        }
    return -1;
}

static void rate(char *buf, size_t n, double secs, long items)
{
    if (secs <= 0 || items <= 0) { snprintf(buf, n, "n/a"); return; }
    double ops = (double)items / secs;
    if (ops >= 1e9)      snprintf(buf, n, "%.2f G", ops / 1e9);
    else if (ops >= 1e6) snprintf(buf, n, "%.2f M", ops / 1e6);
    else if (ops >= 1e3) snprintf(buf, n, "%.2f k", ops / 1e3);
    else                 snprintf(buf, n, "%.1f",  ops);
}

static void ratio(char *buf, size_t n, double base, double ours)
{
    if (base <= 0 || ours <= 0) { snprintf(buf, n, "n/a"); return; }
    snprintf(buf, n, "%.2fx", base / ours);
}

static int bestVariant(int m, int o, double *secs)
{
    int best = -1; double b = 0;
    for (int v = 0; v < NVARIANTS; v++) {
        const GpuCell *c = &g_gpu[v][m][o];
        if (!c->built || c->mismatches || c->kernel_s <= 0) continue;
        if (best < 0 || c->kernel_s < b) { b = c->kernel_s; best = v; }
    }
    if (secs) *secs = b;
    return best;
}

static void writeReport(int items, int reps, double elapsed, int complete)
{
    FILE *f = fopen(g_reportPath, "w");
    if (!f) { perror(g_reportPath); return; }

    fprintf(f, "# MPA-OpenCL benchmark report - %s\n\n", g_env.gpuName);
    if (!complete)
        fprintf(f, "> **Partial report.** The run was interrupted or hit its time budget.\n"
                   "> Rows that never ran are marked `n/a`.\n\n");

    fprintf(f, "## 1. System under test\n\n### GPU\n\n");
    fprintf(f, "| Property | Value |\n|---|---|\n");
    fprintf(f, "| Model | %s |\n", g_env.gpuName);
    fprintf(f, "| Vendor | %s |\n", g_env.gpuVendor);
    fprintf(f, "| Global memory | %.2f GiB |\n", (double)g_env.gmem / 1073741824.0);
    fprintf(f, "| Max single allocation | %.2f GiB |\n", (double)g_env.maxAlloc / 1073741824.0);
    fprintf(f, "| Local memory | %.0f KiB |\n", (double)g_env.lmem / 1024.0);
    fprintf(f, "| Global cache | %.0f KiB |\n", (double)g_env.cache / 1024.0);
    fprintf(f, "| Compute units | %u |\n", g_env.cus);
    fprintf(f, "| Max clock | %u MHz |\n", g_env.clockMHz);
    fprintf(f, "| Max work-group size | %zu |\n", g_env.maxWG);
    fprintf(f, "| OpenCL version | %s |\n", g_env.clVersion);
    fprintf(f, "| Driver | %s |\n\n", g_env.driver);

    fprintf(f, "### Host\n\n| Property | Value |\n|---|---|\n");
    fprintf(f, "| CPU | %s |\n", g_env.cpu);
    fprintf(f, "| Logical cores | %ld |\n", g_env.cores);
    fprintf(f, "| OpenMP threads used | %d |\n", g_env.threads);
    fprintf(f, "| RAM | %.1f GB |\n", g_env.ramGB);
    fprintf(f, "| OS | %s |\n", g_env.os);
    fprintf(f, "| Kernel | %s |\n", g_env.kernel);
    fprintf(f, "| Arch | %s |\n", g_env.arch);
    fprintf(f, "| GMP | %s |\n", gmp_version);
    fprintf(f, "| OpenSSL | %s |\n", OpenSSL_version(OPENSSL_VERSION));
    fprintf(f, "| CGBN | %s |\n\n", g_ncgbn ? "cgbn_results.tsv loaded" : "not measured");

    fprintf(f, "## 2. Method\n\n");
    fprintf(f, "- Base workload %d items, scaled down per operator by its cost weight and by modulus size; the exact count is in every row.\n", items);
    fprintf(f, "- %d timed repetitions, **minimum** reported. Two untimed warm-up launches precede them.\n", reps);
    fprintf(f, "- `kernel` times `clEnqueueNDRangeKernel` + `clFinish` only. `e2e` adds the host->device operand writes and the device->host result read.\n");
    fprintf(f, "- CPU baselines run the identical operands (the generator is reseeded per modulus and operation, so every backend sees the same inputs). Temporaries are preallocated outside the timed region, so the figure is the arithmetic, not marshalling.\n");
    fprintf(f, "- Every GPU cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.\n");
    fprintf(f, "- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.\n");
    fprintf(f, "- Total wall time %.1f s.\n\n", elapsed);

    fprintf(f, "## 3. Correctness\n\n");
    fprintf(f, "| Kernel | Configs run | Passed | Mismatched | Build/launch failed |\n|---|---|---|---|---|\n");
    long gtot = 0, gbad = 0;
    for (int v = 0; v < NVARIANTS; v++) {
        long run = 0, pass = 0, bad = 0, fail = 0;
        for (int m = 0; m < NMODULI; m++)
            for (int o = 0; o < NOPS; o++) {
                const GpuCell *c = &g_gpu[v][m][o];
                if (!c->attempted) continue;
                run++;
                if (!c->built) fail++;
                else if (c->mismatches) bad++;
                else pass++;
            }
        if (!run) continue;
        gtot += run; gbad += bad + fail;
        fprintf(f, "| `%s` (%s) | %ld | %ld | %ld | %ld |\n",
                VARIANTS[v].cl, VARIANTS[v].name, run, pass, bad, fail);
    }
    fprintf(f, "\n**%s** - %ld configurations, %ld problems.\n\n",
            gbad ? "FAILURES PRESENT" : "All configurations correct", gtot, gbad);

    fprintf(f, "## 4. Throughput by modulus\n\n");
    fprintf(f, "Operations per second, higher is better. GPU columns are kernel-only.\n\n");
    for (int m = 0; m < NMODULI; m++) {
        fprintf(f, "### %s (%d-bit)\n\n| Operation | items |", MODULI[m].name, MODULI[m].bits);
        for (int v = 0; v < NVARIANTS; v++) fprintf(f, " %s |", VARIANTS[v].name);
        fprintf(f, " GMP 1T | GMP %dT | OpenSSL %dT | CGBN |\n", g_env.threads, g_env.threads);
        fprintf(f, "|---|---|");
        for (int v = 0; v < NVARIANTS; v++) fprintf(f, "---|");
        fprintf(f, "---|---|---|---|\n");

        for (int o = 0; o < NOPS; o++) {
            const CpuRow *cr = &g_cpu[m][o];
            if (!cr->items) continue;
            fprintf(f, "| %s | %ld |", OPS[o].name, cr->items);
            char b[32];
            for (int v = 0; v < NVARIANTS; v++) {
                const GpuCell *c = &g_gpu[v][m][o];
                if (!c->attempted)      fprintf(f, " - |");
                else if (!c->built)     fprintf(f, " build failed |");
                else if (c->mismatches) fprintf(f, " **WRONG** |");
                else { rate(b, sizeof b, c->kernel_s, c->items); fprintf(f, " %s |", b); }
            }
            rate(b, sizeof b, cr->gmp1, cr->items); fprintf(f, " %s |", b);
            rate(b, sizeof b, cr->gmpN, cr->items); fprintf(f, " %s |", b);
            rate(b, sizeof b, cr->ossl, cr->items); fprintf(f, " %s |", b);
            long ci = 0; double cs = cgbnLookup(MODULI[m].name, OPS[o].name, &ci);
            rate(b, sizeof b, cs, ci); fprintf(f, " %s |\n", b);
        }
        fprintf(f, "\n");
    }

    fprintf(f, "## 5. Speedup of the best kernel over each baseline\n\n");
    fprintf(f, "Ratios above 1.00x mean MPA-OpenCL is faster.\n\n");
    for (int m = 0; m < NMODULI; m++) {
        fprintf(f, "### %s (%d-bit)\n\n", MODULI[m].name, MODULI[m].bits);
        fprintf(f, "| Operation | best kernel | vs GMP 1T | vs GMP %dT | vs OpenSSL | vs CGBN | e2e vs GMP %dT |\n",
                g_env.threads, g_env.threads);
        fprintf(f, "|---|---|---|---|---|---|---|\n");
        for (int o = 0; o < NOPS; o++) {
            const CpuRow *cr = &g_cpu[m][o];
            if (!cr->items) continue;
            double bs = 0; int bv = bestVariant(m, o, &bs);
            if (bv < 0) { fprintf(f, "| %s | none correct | n/a | n/a | n/a | n/a | n/a |\n", OPS[o].name); continue; }
            char r1[32], r2[32], r3[32], r4[32], r5[32];
            long ci = 0; double cs = cgbnLookup(MODULI[m].name, OPS[o].name, &ci);
            double cs_scaled = (cs > 0 && ci > 0) ? cs * (double)cr->items / (double)ci : -1;
            ratio(r1, sizeof r1, cr->gmp1, bs);
            ratio(r2, sizeof r2, cr->gmpN, bs);
            ratio(r3, sizeof r3, cr->ossl, bs);
            ratio(r4, sizeof r4, cs_scaled, bs);
            ratio(r5, sizeof r5, cr->gmpN, g_gpu[bv][m][o].e2e_s);
            fprintf(f, "| %s | %s | %s | %s | %s | %s | %s |\n",
                    OPS[o].name, VARIANTS[bv].name, r1, r2, r3, r4, r5);
        }
        fprintf(f, "\n");
    }

    if (!g_ncgbn) {
        fprintf(f, "## 6. CGBN\n\n");
        fprintf(f, "Not measured on this run. CGBN is CUDA-only and is not built into this host.\n");
        fprintf(f, "To populate the CGBN columns, produce `cgbn_results.tsv` next to the binary with one\n");
        fprintf(f, "whitespace-separated row per measurement and re-run:\n\n");
        fprintf(f, "```\n# modulus_name  operation_name  items  seconds\nsecp256k1  MODMUL  20000  0.00123\n```\n\n");
        fprintf(f, "`modulus_name` and `operation_name` must match the spellings used in the tables above.\n\n");
    }

    fprintf(f, "## %d. Raw data\n\n", g_ncgbn ? 6 : 7);
    fprintf(f, "Also written to `%s` for analysis.\n\n", g_csvPath);
    fprintf(f, "```csv\nkind,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches\n");
    FILE *c = fopen(g_csvPath, "w");
    if (c) fprintf(c, "kind,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches\n");
    for (int m = 0; m < NMODULI; m++)
        for (int o = 0; o < NOPS; o++) {
            const CpuRow *cr = &g_cpu[m][o];
            if (!cr->items) continue;
            struct { const char *k; double s; } cpu[] = {
                { "gmp-1t", cr->gmp1 }, { "gmp-nt", cr->gmpN }, { "openssl-nt", cr->ossl } };
            for (int i = 0; i < 3; i++) {
                if (cpu[i].s <= 0) continue;
                fprintf(f, "cpu,%s,%s,%d,%s,%ld,%.9f,%.3f,0\n", cpu[i].k, MODULI[m].name,
                        MODULI[m].bits, OPS[o].name, cr->items, cpu[i].s, (double)cr->items/cpu[i].s);
                if (c) fprintf(c, "cpu,%s,%s,%d,%s,%ld,%.9f,%.3f,0\n", cpu[i].k, MODULI[m].name,
                        MODULI[m].bits, OPS[o].name, cr->items, cpu[i].s, (double)cr->items/cpu[i].s);
            }
            for (int v = 0; v < NVARIANTS; v++) {
                const GpuCell *g = &g_gpu[v][m][o];
                if (!g->attempted || !g->built) continue;
                fprintf(f, "gpu-kernel,%s,%s,%d,%s,%ld,%.9f,%.3f,%ld\n", VARIANTS[v].name,
                        MODULI[m].name, MODULI[m].bits, OPS[o].name, g->items, g->kernel_s,
                        (double)g->items/g->kernel_s, g->mismatches);
                fprintf(f, "gpu-e2e,%s,%s,%d,%s,%ld,%.9f,%.3f,%ld\n", VARIANTS[v].name,
                        MODULI[m].name, MODULI[m].bits, OPS[o].name, g->items, g->e2e_s,
                        (double)g->items/g->e2e_s, g->mismatches);
                if (c) {
                    fprintf(c, "gpu-kernel,%s,%s,%d,%s,%ld,%.9f,%.3f,%ld\n", VARIANTS[v].name,
                            MODULI[m].name, MODULI[m].bits, OPS[o].name, g->items, g->kernel_s,
                            (double)g->items/g->kernel_s, g->mismatches);
                    fprintf(c, "gpu-e2e,%s,%s,%d,%s,%ld,%.9f,%.3f,%ld\n", VARIANTS[v].name,
                            MODULI[m].name, MODULI[m].bits, OPS[o].name, g->items, g->e2e_s,
                            (double)g->items/g->e2e_s, g->mismatches);
                }
            }
        }
    fprintf(f, "```\n");
    if (c) fclose(c);
    fclose(f);
}
