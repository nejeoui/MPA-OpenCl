#define CL_TARGET_OPENCL_VERSION 120
#define _POSIX_C_SOURCE 200809L
#ifdef __APPLE__
#ifndef _DARWIN_C_SOURCE
#define _DARWIN_C_SOURCE
#endif
#endif
#define OPENSSL_API_COMPAT 0x10100000L
#define OPENSSL_SUPPRESS_DEPRECATED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <gmp.h>
#include <openssl/ec.h>
#include <openssl/ecdsa.h>
#include <openssl/obj_mac.h>
#include <openssl/bn.h>
#include <openssl/rand.h>
#include <openssl/sha.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#include <sys/sysctl.h>
#else
#include <CL/cl.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#define T 8
#define GREEN "\x1b[32m"
#define RED   "\x1b[31m"
#define DIM   "\x1b[2m"
#define OFF   "\x1b[0m"

#define CHECK(e) do { cl_int _e=(e); if(_e!=CL_SUCCESS){ \
  fprintf(stderr,"%s:%d OpenCL error %d\n",__FILE__,__LINE__,(int)_e); exit(2);} } while(0)

static const char *P256_P  = "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF";
static const char *P256_N  = "FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551";
static const char *P256_GX = "6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296";
static const char *P256_GY = "4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5";

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
    int n = 0; size_t sz = sizeof(n);
    if (sysctlbyname("hw.logicalcpu", &n, &sz, NULL, 0) == 0 && n > 0) return n;
    return -1;
#elif defined(_SC_NPROCESSORS_ONLN)
    long n = sysconf(_SC_NPROCESSORS_ONLN); return n > 0 ? n : -1;
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

static const char *devTypeName(cl_device_type t)
{
    if (t & CL_DEVICE_TYPE_GPU) return "GPU";
    if (t & CL_DEVICE_TYPE_ACCELERATOR) return "ACCELERATOR";
    if (t & CL_DEVICE_TYPE_CPU) return "CPU";
    return "OTHER";
}

static void pickDevice(cl_device_id *outDev)
{
    cl_platform_id plats[16];
    cl_uint nplat = 0, nd = 0;
    cl_device_type order[3] = { CL_DEVICE_TYPE_GPU, CL_DEVICE_TYPE_ACCELERATOR,
                                CL_DEVICE_TYPE_CPU };
    const char *want = getenv("MPA_DEVICE_TYPE");
    int t0 = 0, t1 = 3, strict = 0;
    cl_device_id d;

    CHECK(clGetPlatformIDs(16, plats, &nplat));
    if (!nplat) { fprintf(stderr, "no OpenCL platform\n"); exit(2); }
    if (want) {
        strict = 1;
        if      (!strcmp(want, "gpu")) { t0 = 0; t1 = 1; }
        else if (!strcmp(want, "cpu")) { t0 = 2; t1 = 3; }
        else { t0 = 0; t1 = 3; strict = 0; }
    }
    for (int t = t0; t < t1; t++)
        for (cl_uint i = 0; i < nplat; i++)
            if (clGetDeviceIDs(plats[i], order[t], 1, &d, &nd) == CL_SUCCESS && nd > 0) {
                *outDev = d; return;
            }
    fprintf(stderr, strict ? "no %s device found; refusing to fall back\n"
                           : "no OpenCL device found\n", want ? want : "");
    exit(2);
}

static void mpzToWords(const mpz_t z, uint32_t *buf)
{
    mpz_t t; mpz_init_set(t, z);
    for (int i = T - 1; i >= 0; i--) {
        buf[i] = (uint32_t)(mpz_get_ui(t) & 0xFFFFFFFFu);
        mpz_tdiv_q_2exp(t, t, 32);
    }
    mpz_clear(t);
}

static void bnToWords(const BIGNUM *b, uint32_t *buf)
{
    char *h = BN_bn2hex(b);
    mpz_t z; mpz_init_set_str(z, h, 16);
    mpzToWords(z, buf);
    mpz_clear(z); OPENSSL_free(h);
}

int main(int argc, char **argv)
{
    size_t items = 4096;
    int reps = 9, tamperPct = 10;
    const char *kernelFile = "ecdsaKernel_opt.cl";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--items") && i+1 < argc) items = strtoul(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--reps")  && i+1 < argc) reps  = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--tamper")&& i+1 < argc) tamperPct = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--kernel")&& i+1 < argc) kernelFile = argv[++i];
        else if (!strcmp(argv[i], "--baseline")) kernelFile = "ecdsaKernel.cl";
        else { fprintf(stderr, "usage: %s [--items N] [--reps R] [--tamper PCT]\n"
                               "          [--kernel FILE] [--baseline]\n", argv[0]); return 2; }
    }
    if (reps < 1) reps = 1;

    mpz_t P, N, GX, GY, R2P, R2N, base, inv, tmp;
    mpz_inits(P, N, GX, GY, R2P, R2N, base, inv, tmp, NULL);
    mpz_set_str(P, P256_P, 16);
    mpz_set_str(N, P256_N, 16);
    mpz_set_str(GX, P256_GX, 16);
    mpz_set_str(GY, P256_GY, 16);

    mpz_ui_pow_ui(base, 2, 32);
    mpz_invert(inv, P, base); mpz_sub(tmp, base, inv);
    unsigned long mprime_p = mpz_get_ui(tmp);
    mpz_invert(inv, N, base); mpz_sub(tmp, base, inv);
    unsigned long mprime_n = mpz_get_ui(tmp);

    mpz_ui_pow_ui(R2P, 2, 512); mpz_mod(R2P, R2P, P);
    mpz_ui_pow_ui(R2N, 2, 512); mpz_mod(R2N, R2N, N);

    uint32_t par[9 * T];
    memset(par, 0, sizeof(par));
    mpzToWords(P, par);
    mpzToWords(R2P, par + T);
    mpzToWords(N, par + 2 * T);
    mpzToWords(R2N, par + 3 * T);
    { /* G in Montgomery form */
        mpz_t gm; mpz_init(gm);
        mpz_mul_2exp(gm, GX, 256); mpz_mod(gm, gm, P); mpzToWords(gm, par + 4 * T);
        mpz_mul_2exp(gm, GY, 256); mpz_mod(gm, gm, P); mpzToWords(gm, par + 5 * T);
        mpz_clear(gm);
    }
    par[6 * T]     = (uint32_t)mprime_p;
    par[6 * T + 1] = (uint32_t)mprime_n;
    /* r + n can exceed p only when p - n is large enough; P-256 has p > n. */
    {
        mpz_t d; mpz_init(d); mpz_sub(d, P, N);
        par[6 * T + 2] = (mpz_sgn(d) > 0) ? 1u : 0u;
        mpz_clear(d);
    }
    mpzToWords(N, par + 8 * T);

    printf("workload  : %zu P-256 signatures, %d%% deliberately tampered\n",
           items, tamperPct);

    uint32_t *hSig = calloc(items * 2 * T, 4);
    uint32_t *hMsg = calloc(items * T, 4);
    uint32_t *hPub = calloc(items * 2 * T, 4);
    unsigned char *hVal = calloc(items, 1);
    unsigned char *expect = calloc(items, 1);

    EC_KEY **keys = malloc(sizeof(EC_KEY *) * items);
    unsigned char (*digests)[32] = malloc(32 * items);
    ECDSA_SIG **sigs = malloc(sizeof(ECDSA_SIG *) * items);

    double tgen = now_s();
    EC_GROUP *grp = EC_GROUP_new_by_curve_name(NID_X9_62_prime256v1);
    BN_CTX *bctx = BN_CTX_new();
    long nTamper = 0;

    for (size_t i = 0; i < items; i++) {
        EC_KEY *k = EC_KEY_new_by_curve_name(NID_X9_62_prime256v1);
        if (!k || !EC_KEY_generate_key(k)) { fprintf(stderr, "keygen failed\n"); return 2; }
        keys[i] = k;

        RAND_bytes(digests[i], 32);
        ECDSA_SIG *sg = ECDSA_do_sign(digests[i], 32, k);
        if (!sg) { fprintf(stderr, "sign failed\n"); return 2; }

        const BIGNUM *r, *s;
        ECDSA_SIG_get0(sg, &r, &s);

        BIGNUM *rr = BN_dup(r), *ss = BN_dup(s);
        if ((int)(i % 100) < tamperPct) {
            /* flip one bit of r: the signature must now be rejected */
            BN_set_bit(rr, (int)(i % 200));
            BN_clear_bit(rr, (int)((i + 1) % 200));
            nTamper++;
        }
        ECDSA_SIG *sg2 = ECDSA_SIG_new();
        ECDSA_SIG_set0(sg2, BN_dup(rr), BN_dup(ss));
        expect[i] = (unsigned char)(ECDSA_do_verify(digests[i], 32, sg2, k) == 1);
        sigs[i] = sg2;

        bnToWords(rr, hSig + i * 2 * T);
        bnToWords(ss, hSig + i * 2 * T + T);
        BN_free(rr); BN_free(ss);
        ECDSA_SIG_free(sg);

        { /* e = digest as a 256-bit integer */
            mpz_t e; mpz_init(e);
            mpz_import(e, 32, 1, 1, 1, 0, digests[i]);
            mpzToWords(e, hMsg + i * T);
            mpz_clear(e);
        }
        { /* public key affine coordinates */
            BIGNUM *qx = BN_new(), *qy = BN_new();
            EC_POINT_get_affine_coordinates(grp, EC_KEY_get0_public_key(k), qx, qy, bctx);
            bnToWords(qx, hPub + i * 2 * T);
            bnToWords(qy, hPub + i * 2 * T + T);
            BN_free(qx); BN_free(qy);
        }
    }
    tgen = now_s() - tgen;
    printf("generation: %.2f s with OpenSSL (%ld tampered, %ld valid)\n",
           tgen, nTamper, (long)items - nTamper);

    cl_device_id dev;
    pickDevice(&dev);
    char dname[256] = {0}; cl_device_type dtype = 0; cl_uint cus = 0;
    clGetDeviceInfo(dev, CL_DEVICE_NAME, sizeof(dname), dname, NULL);
    clGetDeviceInfo(dev, CL_DEVICE_TYPE, sizeof(dtype), &dtype, NULL);
    clGetDeviceInfo(dev, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cus), &cus, NULL);
    const long cores = hwCoreCount();
    printf("opencl    : [%s] %s (%u compute units)\n", devTypeName(dtype), dname, cus);
    printf("cpu       : %d thread(s), %ld core(s) online\n", nthreads(), cores);
    printf("openssl   : %s\n", OpenSSL_version(OPENSSL_VERSION));
    printf("kernel    : %s\n\n", kernelFile);

    cl_int err;
    cl_context ctx = clCreateContext(NULL, 1, &dev, NULL, NULL, &err); CHECK(err);
    cl_command_queue q = clCreateCommandQueue(ctx, dev, 0, &err); CHECK(err);

    FILE *fp = fopen(kernelFile, "rb");
    if (!fp) { fprintf(stderr, "cannot open %s\n", kernelFile); return 2; }
    char *src = malloc(1 << 20);
    size_t srcLen = fread(src, 1, 1 << 20, fp);
    fclose(fp);

    cl_program prog = clCreateProgramWithSource(ctx, 1, (const char **)&src, &srcLen, &err);
    CHECK(err);
    if (clBuildProgram(prog, 1, &dev, "-DMPA_UNROLL=1", NULL, NULL) != CL_SUCCESS) {
        size_t ln = 0;
        clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &ln);
        char *log = malloc(ln + 1);
        clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, ln, log, NULL);
        log[ln] = 0;
        fprintf(stderr, "kernel build failed:\n%s\n", log);
        return 2;
    }
    cl_kernel kern = clCreateKernel(prog, "ecdsaVerify", &err); CHECK(err);

    cl_mem dSig = clCreateBuffer(ctx, CL_MEM_READ_ONLY, items*2*T*4, NULL, &err); CHECK(err);
    cl_mem dMsg = clCreateBuffer(ctx, CL_MEM_READ_ONLY, items*T*4, NULL, &err); CHECK(err);
    cl_mem dPub = clCreateBuffer(ctx, CL_MEM_READ_ONLY, items*2*T*4, NULL, &err); CHECK(err);
    cl_mem dVal = clCreateBuffer(ctx, CL_MEM_WRITE_ONLY, items, NULL, &err); CHECK(err);
    cl_mem dPar = clCreateBuffer(ctx, CL_MEM_READ_ONLY, sizeof(par), NULL, &err); CHECK(err);

    CHECK(clSetKernelArg(kern, 0, sizeof(cl_mem), &dSig));
    CHECK(clSetKernelArg(kern, 1, sizeof(cl_mem), &dMsg));
    CHECK(clSetKernelArg(kern, 2, sizeof(cl_mem), &dPub));
    CHECK(clSetKernelArg(kern, 3, sizeof(cl_mem), &dVal));
    CHECK(clSetKernelArg(kern, 4, sizeof(cl_mem), &dPar));

    CHECK(clEnqueueWriteBuffer(q, dSig, CL_TRUE, 0, items*2*T*4, hSig, 0, NULL, NULL));
    CHECK(clEnqueueWriteBuffer(q, dMsg, CL_TRUE, 0, items*T*4, hMsg, 0, NULL, NULL));
    CHECK(clEnqueueWriteBuffer(q, dPub, CL_TRUE, 0, items*2*T*4, hPub, 0, NULL, NULL));
    CHECK(clEnqueueWriteBuffer(q, dPar, CL_TRUE, 0, sizeof(par), par, 0, NULL, NULL));

    CHECK(clEnqueueNDRangeKernel(q, kern, 1, NULL, &items, NULL, 0, NULL, NULL));
    CHECK(clFinish(q));
    CHECK(clEnqueueReadBuffer(q, dVal, CL_TRUE, 0, items, hVal, 0, NULL, NULL));

    long falseAccept = 0, falseReject = 0;
    for (size_t i = 0; i < items; i++) {
        if (hVal[i] && !expect[i]) falseAccept++;
        if (!hVal[i] && expect[i]) falseReject++;
    }
    printf("agreement with OpenSSL: ");
    if (falseAccept || falseReject)
        printf(RED "%ld false accepts, %ld false rejects" OFF "\n", falseAccept, falseReject);
    else
        printf(GREEN "exact on all %zu signatures" OFF "\n", items);
    if (falseAccept)
        printf(RED "  FALSE ACCEPTS ARE A SECURITY FAILURE - results below are void\n" OFF);
    printf("\n");

    for (int w = 0; w < 2; w++) {
        CHECK(clEnqueueNDRangeKernel(q, kern, 1, NULL, &items, NULL, 0, NULL, NULL));
        CHECK(clFinish(q));
    }
    double *tk = malloc(sizeof(double) * reps);
    double *te = malloc(sizeof(double) * reps);
    for (int r = 0; r < reps; r++) {
        double t0 = now_s();
        CHECK(clEnqueueNDRangeKernel(q, kern, 1, NULL, &items, NULL, 0, NULL, NULL));
        CHECK(clFinish(q));
        tk[r] = now_s() - t0;

        double t1 = now_s();
        CHECK(clEnqueueWriteBuffer(q, dSig, CL_FALSE, 0, items*2*T*4, hSig, 0, NULL, NULL));
        CHECK(clEnqueueWriteBuffer(q, dMsg, CL_FALSE, 0, items*T*4, hMsg, 0, NULL, NULL));
        CHECK(clEnqueueWriteBuffer(q, dPub, CL_FALSE, 0, items*2*T*4, hPub, 0, NULL, NULL));
        CHECK(clEnqueueNDRangeKernel(q, kern, 1, NULL, &items, NULL, 0, NULL, NULL));
        CHECK(clEnqueueReadBuffer(q, dVal, CL_TRUE, 0, items, hVal, 0, NULL, NULL));
        CHECK(clFinish(q));
        te[r] = now_s() - t1;
    }
    double clK = minimum(tk, reps), clE = minimum(te, reps);

    const int NT = nthreads();
    double *to = malloc(sizeof(double) * reps);
    for (int r = 0; r < reps; r++) {
        double t0 = now_s();
        for (size_t i = 0; i < items; i++)
            ECDSA_do_verify(digests[i], 32, sigs[i], keys[i]);
        to[r] = now_s() - t0;
    }
    double o1 = minimum(to, reps);

    double oN = o1;
#ifdef _OPENMP
    if (NT > 1) {
        for (int r = 0; r < reps; r++) {
            double t0 = now_s();
#pragma omp parallel for schedule(static)
            for (size_t i = 0; i < items; i++)
                ECDSA_do_verify(digests[i], 32, sigs[i], keys[i]);
            to[r] = now_s() - t0;
        }
        oN = minimum(to, reps);
    }
#endif

    printf("  %-26s %11s %14s %11s\n", "backend", "time ms", "verif/s", "vs OpenSSL");
    printf("  ------------------------------------------------------------------\n");
    double best = o1 < oN ? o1 : oN;
    struct { const char *n; double t; } rows[] = {
        { "OpenSSL, 1 thread",   o1  },
        { "OpenSSL, all threads", oN },
        { "OpenCL kernel only",  clK },
        { "OpenCL end-to-end",   clE },
    };
    for (int i = 0; i < 4; i++)
        printf("  %-26s %11.3f %14.3e %10.2fx\n", rows[i].n, rows[i].t * 1e3,
               (double)items / rows[i].t, best / rows[i].t);

    printf("\n  latency for one batch of %zu: %.2f ms end-to-end\n", items, clE * 1e3);
    printf("  per-signature amortised    : %.2f us\n", clE * 1e6 / (double)items);
    if (nthreads() == 1 && cores > 1)
        printf(RED "  WARNING: OpenSSL rows are single-threaded of %ld cores\n" OFF, cores);
    if (!(dtype & CL_DEVICE_TYPE_GPU))
        printf(RED "  NOTE: OpenCL device is not a GPU\n" OFF);

    return (falseAccept || falseReject) ? 1 : 0;
}
