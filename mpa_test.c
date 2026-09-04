#include "mpa_ref.h"

typedef struct { long checked, failed; } Tally;

static void runBatch(cl_context ctx, cl_command_queue q, cl_kernel kern,
                     cl_device_id dev,
                     const Op *op, const Modulus *mod, const Variant *var,
                     int T, int items, const mpz_t p, unsigned long m_prime,
                     Tally *tally)
{
    const int    wbits    = var->wbits;
    const size_t wsz      = (size_t)(wbits / 8);
    const int    bits     = T * wbits;
    const int    outWords = op->wide ? 2 * T : T;

    void *hA   = calloc((size_t)items * T,        wsz);
    void *hB   = calloc((size_t)items * T,        wsz);
    void *hC   = calloc((size_t)items * outWords, wsz);
    void *hExp = calloc((size_t)items * outWords, wsz);
    void *hP   = calloc((size_t)T * 2,            wsz);
    cl_int opBuf[4];

    mpz_t a, b, e;
    mpz_inits(a, b, e, NULL);

    for (int j = 0; j < items; j++) {
        makeCase(j, op->op, op->modular, p, bits, mod, a, b);
        mpzToWords(a, hA, (size_t)j * T, T, wbits);
        mpzToWords(b, hB, (size_t)j * T, T, wbits);
        reference(op->op, a, b, p, bits, e);
        mpzToWords(e, hExp, (size_t)j * outWords, outWords, wbits);
    }
    mpzToWords(p, hP, 0, T, wbits);
    {
        mpz_t r2;
        mpz_init(r2);
        computeR2(r2, p, bits);
        mpzToWords(r2, hP, (size_t)T, T, wbits);
        mpz_clear(r2);
    }

    opBuf[0] = op->op;
    opBuf[1] = wbits;
    opBuf[2] = bits;
    opBuf[3] = (cl_int)(uint32_t)m_prime;

    cl_int err;
    cl_mem dA = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  (size_t)items*T*wsz, NULL, &err);        CHECK(err);
    cl_mem dB = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  (size_t)items*T*wsz, NULL, &err);        CHECK(err);
    cl_mem dC = clCreateBuffer(ctx, CL_MEM_READ_WRITE, (size_t)items*outWords*wsz, NULL, &err); CHECK(err);
    cl_mem dO = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  sizeof(opBuf), NULL, &err);              CHECK(err);
    cl_mem dP = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  (size_t)T*2*wsz, NULL, &err);              CHECK(err);

    CHECK(clEnqueueWriteBuffer(q, dA, CL_TRUE, 0, (size_t)items*T*wsz, hA, 0, NULL, NULL));
    CHECK(clEnqueueWriteBuffer(q, dB, CL_TRUE, 0, (size_t)items*T*wsz, hB, 0, NULL, NULL));
    CHECK(clEnqueueWriteBuffer(q, dO, CL_TRUE, 0, sizeof(opBuf), opBuf, 0, NULL, NULL));
    CHECK(clEnqueueWriteBuffer(q, dP, CL_TRUE, 0, (size_t)T*2*wsz, hP, 0, NULL, NULL));

    CHECK(clSetKernelArg(kern, 0, sizeof(cl_mem), &dA));
    CHECK(clSetKernelArg(kern, 1, sizeof(cl_mem), &dB));
    CHECK(clSetKernelArg(kern, 2, sizeof(cl_mem), &dC));
    CHECK(clSetKernelArg(kern, 3, sizeof(cl_mem), &dO));
    CHECK(clSetKernelArg(kern, 4, sizeof(cl_mem), &dP));

    size_t global = (size_t)items;
    CHECK(clEnqueueNDRangeKernel(q, kern, 1, NULL, &global, NULL, 0, NULL, NULL));
    CHECK(clFinish(q));
    CHECK(clEnqueueReadBuffer(q, dC, CL_TRUE, 0, (size_t)items*outWords*wsz, hC, 0, NULL, NULL));

    long bad = 0;
    int  firstBad = -1;
    for (int j = 0; j < items; j++) {
        int ok = 1;
        for (int i = 0; i < outWords; i++) {
            if (loadWord(hC, (size_t)j*outWords + i, wbits) !=
                loadWord(hExp, (size_t)j*outWords + i, wbits)) { ok = 0; break; }
        }
        if (!ok) { bad++; if (firstBad < 0) firstBad = j; }
    }

    tally->checked += items;
    tally->failed  += bad;

    printf("  %-26s %-16s %-8s T=%-3d %8d cases  %s%s" OFF "\n",
           op->name, mod->name, var->name, T, items,
           bad ? RED : GREEN, bad ? "FAIL" : "ok");

    if (bad) {
        printf("      %ld/%d mismatched (first at case %d%s)\n",
               bad, items, firstBad, firstBad < EDGE_CASES ? ", a directed edge case" : "");
        if (g_verbose) {
            mpz_t got, want, ga, gb;
            mpz_inits(got, want, ga, gb, NULL);
            wordsToMpz(ga,   hA,   (size_t)firstBad*T,        T,        wbits);
            wordsToMpz(gb,   hB,   (size_t)firstBad*T,        T,        wbits);
            wordsToMpz(got,  hC,   (size_t)firstBad*outWords, outWords, wbits);
            wordsToMpz(want, hExp, (size_t)firstBad*outWords, outWords, wbits);
            gmp_printf("      a    = %Zx\n      b    = %Zx\n"
                       "      got  = %Zx\n      want = %Zx\n", ga, gb, got, want);
            mpz_clears(got, want, ga, gb, NULL);
        }
    }

    clReleaseMemObject(dA); clReleaseMemObject(dB); clReleaseMemObject(dC);
    clReleaseMemObject(dO); clReleaseMemObject(dP);
    mpz_clears(a, b, e, NULL);
    free(hA); free(hB); free(hC); free(hExp); free(hP);
}

int main(int argc, char **argv)
{
    int items = 2000;
    int onlyWidth = 0;
    unsigned long seed = 12345;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--items")   && i+1 < argc) items = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed") && i+1 < argc) seed = strtoul(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--width") && i+1 < argc) onlyWidth = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--verbose")) g_verbose = 1;
        else { fprintf(stderr, "usage: %s [--items N] [--seed S] [--width 8|16|32] [--verbose]\n", argv[0]); return 2; }
    }
    if (items < EDGE_CASES) items = EDGE_CASES;
    rngState = seed;

    cl_platform_id plat;
    cl_device_id   dev;
    pickDevice(&plat, &dev);

    char nameBuf[256] = {0}, verBuf[128] = {0};
    cl_device_type dtype = 0;
    clGetDeviceInfo(dev, CL_DEVICE_NAME, sizeof(nameBuf), nameBuf, NULL);
    clGetDeviceInfo(dev, CL_DEVICE_VERSION, sizeof(verBuf), verBuf, NULL);
    clGetDeviceInfo(dev, CL_DEVICE_TYPE, sizeof(dtype), &dtype, NULL);
    printf("device : [%s] %s\n", devTypeName(dtype), nameBuf);
    printf("version: %s\n", verBuf);
    printf("items  : %d per (op, modulus, width), seed %lu\n\n", items, seed);

    cl_int err;
    cl_context ctx = clCreateContext(NULL, 1, &dev, NULL, NULL, &err); CHECK(err);
    cl_command_queue q = clCreateCommandQueue(ctx, dev, 0, &err);      CHECK(err);

    Tally tally = { 0, 0 };
    int programs = 0;

    for (int v = 0; v < NVARIANTS; v++) {
        const Variant *var = &VARIANTS[v];
        if (onlyWidth && var->wbits != onlyWidth) continue;

        size_t srcLen;
        char  *src = readFile(var->cl, &srcLen);

        for (int m = 0; m < NMODULI; m++) {
            const Modulus *mod = &MODULI[m];
            int T = mod->bits / var->wbits;

            mpz_t p, base, inv, mp;
            mpz_inits(p, base, inv, mp, NULL);
            if (mpz_set_str(p, mod->hex, 16) != 0) {
                fprintf(stderr, "bad modulus literal for %s\n", mod->name); return 2;
            }
            if ((int)mpz_sizeinbase(p, 2) != mod->bits || mpz_even_p(p)) {
                fprintf(stderr, "%s: literal is %d bits (expected %d) or is even\n",
                        mod->name, (int)mpz_sizeinbase(p, 2), mod->bits);
                return 2;
            }

            mpz_ui_pow_ui(base, 2, (unsigned long)var->wbits);
            if (mpz_invert(inv, p, base) == 0) {
                fprintf(stderr, "%s is even; no Montgomery inverse\n", mod->name); return 2;
            }
            mpz_sub(mp, base, inv);
            unsigned long m_prime = mpz_get_ui(mp);

            char opts[256];
            snprintf(opts, sizeof(opts), "-I. -DWORDLENGTH_T=%d %s", T, var->flags);

            cl_program prog = clCreateProgramWithSource(ctx, 1, (const char **)&src, &srcLen, &err);
            CHECK(err);
            cl_int berr = clBuildProgram(prog, 1, &dev, opts, NULL, NULL);
            if (berr != CL_SUCCESS) {
                size_t logn = 0;
                clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &logn);
                char *log = malloc(logn + 1);
                clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, logn, log, NULL);
                log[logn] = 0;
                fprintf(stderr, "build failed for %s (%s):\n%s\n", var->cl, opts, log);
                free(log);
                return 2;
            }
            cl_kernel kern = clCreateKernel(prog, "mpaKernel", &err); CHECK(err);
            programs++;

            printf(DIM "-- %s / %s / %d-bit words / m'=0x%lx --" OFF "\n",
                   mod->name, var->name, var->wbits, m_prime);

            for (int o = 0; o < NOPS; o++) {
                if (OPS[o].ext && !var->ext) continue;
                int n = items;
                int div = OPS[o].cost * (T > 8 ? T / 8 : 1);
                if (div > 1) n = items / div;
                if (n < EDGE_CASES + 16) n = EDGE_CASES + 16;
                runBatch(ctx, q, kern, dev, &OPS[o], mod, var, T, n, p, m_prime, &tally);
            }

            printf("\n");
            clReleaseKernel(kern);
            clReleaseProgram(prog);
            mpz_clears(p, base, inv, mp, NULL);
        }
        free(src);
    }

    clReleaseCommandQueue(q);
    clReleaseContext(ctx);

    printf("================================================================\n");
    printf("%d kernel builds, %ld cases checked, %s%ld failed" OFF "\n",
           programs, tally.checked, tally.failed ? RED : GREEN, tally.failed);
    printf("================================================================\n");
    return tally.failed ? 1 : 0;
}
