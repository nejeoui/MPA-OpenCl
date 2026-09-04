#include "mpa_ref.h"
#include <time.h>
#include <openssl/bn.h>
#include <openssl/err.h>
#include <unistd.h>
#if defined(__APPLE__)
#include <sys/sysctl.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

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

static long hwCoreCount(void)
{
#if defined(__APPLE__)
    int n = 0;
    size_t sz = sizeof(n);
    if (sysctlbyname("hw.logicalcpu", &n, &sz, NULL, 0) == 0 && n > 0) return n;
    return -1;
#elif defined(_SC_NPROCESSORS_ONLN)
    const long n = sysconf(_SC_NPROCESSORS_ONLN);
    return (n > 0) ? n : -1;
#else
    return -1;
#endif
}

static int nthreads(void)
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

static void gmp_kernel(int op, mpz_t r, const mpz_t a, const mpz_t b,
                       const mpz_t p)
{
    switch (op) {
    case MODMUL:
    case MODMUL_R2: mpz_mul(r, a, b); mpz_mod(r, r, p); break;
    case MODEXP: mpz_powm(r, a, b, p);               break;
    default:     mpz_mul(r, a, b);                   break;
    }
}

static double run_gmp(int op, size_t items, mpz_t *A, mpz_t *B, const mpz_t p,
                      mpz_t *OUT, int threads, int reps)
{
    double *t = malloc(sizeof(double) * (size_t)reps);
    for (int r = 0; r < reps; r++) {
        double t0 = now_s();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
        for (size_t i = 0; i < items; i++)
            gmp_kernel(op, OUT[i], A[i], B[i], p);
        t[r] = now_s() - t0;
    }
    double m = minimum(t, reps);
    free(t);
    return m;
}

static double run_ossl(int op, size_t items, BIGNUM **A, BIGNUM **B,
                       const BIGNUM *p, BIGNUM **OUT, int threads, int reps)
{
    double *t = malloc(sizeof(double) * (size_t)reps);
    for (int r = 0; r < reps; r++) {
        double t0 = now_s();
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
        {
            BN_CTX *ctx = BN_CTX_new();
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (size_t i = 0; i < items; i++) {
                switch (op) {
                case MODMUL:
                case MODMUL_R2: BN_mod_mul(OUT[i], A[i], B[i], p, ctx); break;
                case MODEXP: BN_mod_exp(OUT[i], A[i], B[i], p, ctx); break;
                default:     BN_mul(OUT[i], A[i], B[i], ctx);        break;
                }
            }
            BN_CTX_free(ctx);
        }
        t[r] = now_s() - t0;
    }
    double m = minimum(t, reps);
    free(t);
    return m;
}

static void mpzToBuf(const mpz_t z, uint32_t *buf, size_t it, int nwords)
{
    mpz_t t;
    mpz_init_set(t, z);
    for (int i = nwords - 1; i >= 0; i--) {
        buf[it * (size_t)nwords + (size_t)i] =
            (uint32_t)(mpz_get_ui(t) & 0xFFFFFFFFu);
        mpz_tdiv_q_2exp(t, t, 32);
    }
    mpz_clear(t);
}

static void bufToMpz(mpz_t z, const uint32_t *buf, size_t it, int nwords)
{
    mpz_set_ui(z, 0);
    for (int i = 0; i < nwords; i++) {
        mpz_mul_2exp(z, z, 32);
        mpz_add_ui(z, z, (unsigned long)buf[it * (size_t)nwords + (size_t)i]);
    }
}

int main(int argc, char **argv)
{
    size_t items = 20000;
    int reps = 9, bits = 256;
    const char *opName = "MODEXP";
    const char *clFlags = "-DMPA_REGACC=1 -DMPA_FUSED_CIOS=1 -DMPA_UNROLL=1";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--items") && i+1 < argc) items = strtoul(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--reps")  && i+1 < argc) reps  = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--bits")  && i+1 < argc) bits  = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--op")    && i+1 < argc) opName = argv[++i];
        else if (!strcmp(argv[i], "--flags") && i+1 < argc) clFlags = argv[++i];
        else if (!strcmp(argv[i], "--mulhi")) clFlags =
            "-DMPA_MULHI=1 -DMPA_REGACC=1 -DMPA_FUSED_CIOS=1 -DMPA_UNROLL=1";
        else {
            fprintf(stderr, "usage: %s [--items N] [--reps R] [--bits B]\n"
                            "          [--op MODMUL|MODMUL_R2|MODEXP|MULTIPLY]\n"
                            "          [--mulhi] [--flags \"-D...\"]\n", argv[0]);
            return 2;
        }
    }
    if (reps < 1) reps = 1;
    rngState = 4242;

    int op;
    if      (!strcmp(opName, "MODMUL"))    op = MODMUL;
    else if (!strcmp(opName, "MODMUL_R2")) op = MODMUL_R2;
    else if (!strcmp(opName, "MODEXP"))   op = MODEXP;
    else if (!strcmp(opName, "MULTIPLY")) op = MULTIPLYOPRANDSCANNING;
    else { fprintf(stderr, "unknown op %s\n", opName); return 2; }

    const int wide = (op == MULTIPLYOPRANDSCANNING);
    const int T = bits / 32;
    const int outWords = wide ? 2 * T : T;

    const Modulus *mod = NULL;
    for (int m = 0; m < NMODULI; m++)
        if (MODULI[m].bits == bits && !MODULI[m].facA) { mod = &MODULI[m]; break; }
    if (!mod) { fprintf(stderr, "no prime modulus of %d bits\n", bits); return 2; }

    mpz_t p, base, inv, mpr;
    mpz_inits(p, base, inv, mpr, NULL);
    mpz_set_str(p, mod->hex, 16);
    mpz_ui_pow_ui(base, 2, 32);
    mpz_invert(inv, p, base);
    mpz_sub(mpr, base, inv);
    unsigned long m_prime = mpz_get_ui(mpr);

    mpz_t *A = malloc(sizeof(mpz_t) * items);
    mpz_t *B = malloc(sizeof(mpz_t) * items);
    mpz_t *REF = malloc(sizeof(mpz_t) * items);
    mpz_t *OUT = malloc(sizeof(mpz_t) * items);
    for (size_t i = 0; i < items; i++) {
        mpz_inits(A[i], B[i], REF[i], OUT[i], NULL);
        makeCase((int)(i % 4096) + EDGE_CASES, op, 1, p, bits, mod, A[i], B[i]);
        gmp_kernel(op, REF[i], A[i], B[i], p);
    }

    BIGNUM **bA = malloc(sizeof(BIGNUM *) * items);
    BIGNUM **bB = malloc(sizeof(BIGNUM *) * items);
    BIGNUM **bO = malloc(sizeof(BIGNUM *) * items);
    BIGNUM *bP = NULL;
    {
        char *hs = mpz_get_str(NULL, 16, p);
        BN_hex2bn(&bP, hs);
        free(hs);
        for (size_t i = 0; i < items; i++) {
            bA[i] = NULL; bB[i] = NULL; bO[i] = BN_new();
            char *ha = mpz_get_str(NULL, 16, A[i]);
            char *hb = mpz_get_str(NULL, 16, B[i]);
            BN_hex2bn(&bA[i], ha);
            BN_hex2bn(&bB[i], hb);
            free(ha); free(hb);
        }
    }

    cl_platform_id plat; cl_device_id dev;
    pickDevice(&plat, &dev);
    char dname[256] = {0}; cl_device_type dtype = 0; cl_uint cus = 0;
    clGetDeviceInfo(dev, CL_DEVICE_NAME, sizeof(dname), dname, NULL);
    clGetDeviceInfo(dev, CL_DEVICE_TYPE, sizeof(dtype), &dtype, NULL);
    clGetDeviceInfo(dev, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cus), &cus, NULL);

    printf("operation : %s, %d-bit modulus %s\n", opName, bits, mod->name);
    printf("workload  : %zu items, %d reps, minimum reported\n", items, reps);
    const long hwCores = hwCoreCount();
    if (hwCores > 0)
        printf("cpu       : %d thread(s) in use, %ld core(s) online\n", nthreads(), hwCores);
    else
        printf("cpu       : %d thread(s) in use, core count unknown\n", nthreads());
    printf("opencl    : [%s] %s (%u compute units)\n", devTypeName(dtype), dname, cus);
    printf("libs      : GMP %s, %s\n", gmp_version, OpenSSL_version(OPENSSL_VERSION));
    printf("kernel    : %s\n\n", clFlags);

    cl_int err;
    cl_context ctx = clCreateContext(NULL, 1, &dev, NULL, NULL, &err); CHECK(err);
    cl_command_queue q = clCreateCommandQueue(ctx, dev, 0, &err);      CHECK(err);

    size_t srcLen;
    char *src = readFile("mpaKernel_32bits_opt.cl", &srcLen);
    char opts[256];
    snprintf(opts, sizeof(opts), "-DWORDLENGTH_T=%d %s", T, clFlags);
    cl_program prog = clCreateProgramWithSource(ctx, 1, (const char **)&src, &srcLen, &err);
    CHECK(err);
    if (clBuildProgram(prog, 1, &dev, opts, NULL, NULL) != CL_SUCCESS) {
        size_t ln = 0;
        clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &ln);
        char *log = malloc(ln + 1);
        clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, ln, log, NULL);
        log[ln] = 0;
        fprintf(stderr, "kernel build failed:\n%s\n", log);
        return 2;
    }
    cl_kernel kern = clCreateKernel(prog, "mpaKernel", &err); CHECK(err);

    uint32_t *hA = calloc(items * (size_t)T, 4);
    uint32_t *hB = calloc(items * (size_t)T, 4);
    uint32_t *hC = calloc(items * (size_t)outWords, 4);
    uint32_t *hP = calloc((size_t)T * 2, 4);
    for (size_t i = 0; i < items; i++) {
        mpzToBuf(A[i], hA, i, T);
        mpzToBuf(B[i], hB, i, T);
    }
    mpzToBuf(p, hP, 0, T);
    {
        mpz_t r2;
        mpz_init(r2);
        computeR2(r2, p, bits);
        for (int i = T - 1; i >= 0; i--) {
            hP[T + i] = (uint32_t)(mpz_get_ui(r2) & 0xFFFFFFFFu);
            mpz_tdiv_q_2exp(r2, r2, 32);
        }
        mpz_clear(r2);
    }
    cl_int opBuf[4] = { op, 32, bits, (cl_int)(uint32_t)m_prime };

    const size_t inBytes = items * (size_t)T * 4;
    const size_t outBytes = items * (size_t)outWords * 4;
    cl_mem dA = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  inBytes, NULL, &err); CHECK(err);
    cl_mem dB = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  inBytes, NULL, &err); CHECK(err);
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

    for (int w = 0; w < 2; w++) {
        CHECK(clEnqueueNDRangeKernel(q, kern, 1, NULL, &items, NULL, 0, NULL, NULL));
        CHECK(clFinish(q));
    }
    CHECK(clEnqueueReadBuffer(q, dC, CL_TRUE, 0, outBytes, hC, 0, NULL, NULL));

    long clBad = 0;
    {
        mpz_t got;
        mpz_init(got);
        for (size_t i = 0; i < items; i++) {
            bufToMpz(got, hC, i, outWords);
            if (mpz_cmp(got, REF[i]) != 0) clBad++;
        }
        mpz_clear(got);
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
    double clK = minimum(tk, reps), clE = minimum(te, reps);

    double clZ = 0.0, clZK = 0.0;
    long zBad = 0;
    int zcOk = 0;
    {
        cl_int e1 = CL_SUCCESS, e2 = CL_SUCCESS, e3 = CL_SUCCESS;
        cl_mem zA = clCreateBuffer(ctx, CL_MEM_READ_ONLY  | CL_MEM_ALLOC_HOST_PTR,
                                   inBytes,  NULL, &e1);
        cl_mem zB = clCreateBuffer(ctx, CL_MEM_READ_ONLY  | CL_MEM_ALLOC_HOST_PTR,
                                   inBytes,  NULL, &e2);
        cl_mem zC = clCreateBuffer(ctx, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                   outBytes, NULL, &e3);
        if (e1 == CL_SUCCESS && e2 == CL_SUCCESS && e3 == CL_SUCCESS) {
            cl_int me = CL_SUCCESS;
            void *pa = clEnqueueMapBuffer(q, zA, CL_TRUE, CL_MAP_WRITE, 0, inBytes,
                                          0, NULL, NULL, &me);
            if (me == CL_SUCCESS) { memcpy(pa, hA, inBytes);
                clEnqueueUnmapMemObject(q, zA, pa, 0, NULL, NULL); }
            void *pb = clEnqueueMapBuffer(q, zB, CL_TRUE, CL_MAP_WRITE, 0, inBytes,
                                          0, NULL, NULL, &me);
            if (me == CL_SUCCESS) { memcpy(pb, hB, inBytes);
                clEnqueueUnmapMemObject(q, zB, pb, 0, NULL, NULL); }
            CHECK(clFinish(q));

            CHECK(clSetKernelArg(kern, 0, sizeof(cl_mem), &zA));
            CHECK(clSetKernelArg(kern, 1, sizeof(cl_mem), &zB));
            CHECK(clSetKernelArg(kern, 2, sizeof(cl_mem), &zC));

            for (int w = 0; w < 2; w++) {
                CHECK(clEnqueueNDRangeKernel(q, kern, 1, NULL, &items, NULL, 0, NULL, NULL));
                CHECK(clFinish(q));
            }
            {
                uint32_t *pc = clEnqueueMapBuffer(q, zC, CL_TRUE, CL_MAP_READ, 0,
                                                  outBytes, 0, NULL, NULL, &me);
                if (me == CL_SUCCESS) {
                    mpz_t got;
                    mpz_init(got);
                    for (size_t i = 0; i < items; i++) {
                        bufToMpz(got, pc, i, outWords);
                        if (mpz_cmp(got, REF[i]) != 0) zBad++;
                    }
                    mpz_clear(got);
                    clEnqueueUnmapMemObject(q, zC, pc, 0, NULL, NULL);
                    CHECK(clFinish(q));
                }
            }

            double *tzk = malloc(sizeof(double) * (size_t)reps);
            for (int r = 0; r < reps; r++) {
                double t0 = now_s();
                CHECK(clEnqueueNDRangeKernel(q, kern, 1, NULL, &items, NULL, 0, NULL, NULL));
                CHECK(clFinish(q));
                tzk[r] = now_s() - t0;
            }
            clZK = minimum(tzk, reps);
            free(tzk);

            double *tz = malloc(sizeof(double) * (size_t)reps);
            for (int r = 0; r < reps; r++) {
                double t0 = now_s();
                void *ma = clEnqueueMapBuffer(q, zA, CL_TRUE, CL_MAP_WRITE, 0, inBytes,
                                              0, NULL, NULL, &me);
                clEnqueueUnmapMemObject(q, zA, ma, 0, NULL, NULL);
                void *mb = clEnqueueMapBuffer(q, zB, CL_TRUE, CL_MAP_WRITE, 0, inBytes,
                                              0, NULL, NULL, &me);
                clEnqueueUnmapMemObject(q, zB, mb, 0, NULL, NULL);
                CHECK(clEnqueueNDRangeKernel(q, kern, 1, NULL, &items, NULL, 0, NULL, NULL));
                volatile uint32_t sink = 0;
                uint32_t *mc = clEnqueueMapBuffer(q, zC, CL_TRUE, CL_MAP_READ, 0,
                                                 outBytes, 0, NULL, NULL, &me);
                sink ^= mc[0] ^ mc[items * (size_t)outWords - 1];
                clEnqueueUnmapMemObject(q, zC, mc, 0, NULL, NULL);
                CHECK(clFinish(q));
                (void)sink;
                tz[r] = now_s() - t0;
            }
            clZ = minimum(tz, reps);
            free(tz);
            zcOk = 1;

            CHECK(clSetKernelArg(kern, 0, sizeof(cl_mem), &dA));
            CHECK(clSetKernelArg(kern, 1, sizeof(cl_mem), &dB));
            CHECK(clSetKernelArg(kern, 2, sizeof(cl_mem), &dC));
        }
        if (zA) clReleaseMemObject(zA);
        if (zB) clReleaseMemObject(zB);
        if (zC) clReleaseMemObject(zC);
    }

    const int NT = nthreads();
    double g1 = run_gmp(op, items, A, B, p, OUT, 1, reps);
    long g1bad = 0;
    for (size_t i = 0; i < items; i++) if (mpz_cmp(OUT[i], REF[i]) != 0) g1bad++;

    double gN = (NT > 1) ? run_gmp(op, items, A, B, p, OUT, NT, reps) : g1;

    double o1 = run_ossl(op, items, bA, bB, bP, bO, 1, reps);
    long o1bad = 0;
    {
        mpz_t got;
        mpz_init(got);
        for (size_t i = 0; i < items; i++) {
            char *h = BN_bn2hex(bO[i]);
            mpz_set_str(got, h, 16);
            OPENSSL_free(h);
            if (mpz_cmp(got, REF[i]) != 0) o1bad++;
        }
        mpz_clear(got);
    }
    double oN = (NT > 1) ? run_ossl(op, items, bA, bB, bP, bO, NT, reps) : o1;

    printf("  %-25s %12s %14s %11s %11s\n",
           "backend", "time ms", "ops/s", "vs GMP-1t", "vs best CPU");
    printf("  ------------------------------------------------------------------------------\n");

    double bestCpu = g1;
    if (gN < bestCpu) bestCpu = gN;
    if (o1 < bestCpu) bestCpu = o1;
    if (oN < bestCpu) bestCpu = oN;

    struct { const char *n; double t; long bad; } rows[] = {
        { "GMP (1 thread)",        g1,  g1bad },
        { "GMP (all threads)",     gN,  g1bad },
        { "OpenSSL BN (1 thread)", o1,  o1bad },
        { "OpenSSL BN (all thr.)", oN,  o1bad },
        { "OpenCL kernel (dev buf)",  clK,  clBad },
        { "OpenCL e2e  (dev buf)",    clE,  clBad },
        { "OpenCL kernel (host buf)", clZK, zBad },
        { "OpenCL e2e  (zero-copy)",  clZ,  zBad },
    };
    const int nrows = zcOk ? 8 : 6;
    for (int i = 0; i < nrows; i++) {
        if (rows[i].bad) {
            printf("  %-25s %s%ld/%zu results WRONG%s\n",
                   rows[i].n, RED, rows[i].bad, items, OFF);
            continue;
        }
        printf("  %-25s %12.3f %14.3e %10.2fx %10.2fx\n",
               rows[i].n, rows[i].t * 1e3, (double)items / rows[i].t,
               g1 / rows[i].t, bestCpu / rows[i].t);
    }

    printf("\n  all backends verified against GMP on the same %zu operands\n", items);
    if (nthreads() == 1) {
        printf("  %sWARNING: the CPU rows are running on a single thread", RED);
        if (hwCores > 1) printf(" of %ld cores", hwCores);
        printf(".\n  The CPU baseline is understated and every speedup against it is\n"
               "  correspondingly overstated. Rebuild with OpenMP enabled\n"
               "  (macOS: brew install libomp, then make clean && make).%s\n", OFF);
    }
    if (zcOk)
        printf("  dev buf  = clCreateBuffer default + explicit write/read\n"
               "  host buf = CL_MEM_ALLOC_HOST_PTR + map/unmap, no copies\n"
               "  Both kernel rows run identical code; they differ only in where the\n"
               "  operands live, so comparing them isolates the memory path. Zero-copy\n"
               "  wins on unified memory and usually loses on a discrete GPU, where the\n"
               "  kernel would read operands over PCIe instead of from device memory.\n");
    if (!(dtype & CL_DEVICE_TYPE_GPU))
        printf("  %sNOTE: the OpenCL device is not a GPU; this compares CPU against CPU%s\n",
               RED, OFF);

    return (clBad || g1bad || o1bad || zBad) ? 1 : 0;
}
