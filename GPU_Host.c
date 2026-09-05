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
#define MAXDEV 8

typedef struct {
    cl_device_id id;
    cl_device_type type;
    int    isGpu;
    char   name[256], vendor[160], clver[160], driver[160];
    cl_ulong gmem, lmem, cache, maxAlloc;
    cl_uint  cus, clockMHz;
    size_t   maxWG;
} DevInfo;

static DevInfo g_devs[MAXDEV];
static int     g_ndev = 0;
static int     g_primary = 0;
static GpuCell g_gpu[MAXDEV][NVARIANTS][NMODULI][NOPS];
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

static void collectDevice(DevInfo *d)
{
    cl_device_id x = d->id;
    clGetDeviceInfo(x, CL_DEVICE_TYPE,                 sizeof d->type,     &d->type,     NULL);
    clGetDeviceInfo(x, CL_DEVICE_NAME,                 sizeof d->name,      d->name,     NULL);
    clGetDeviceInfo(x, CL_DEVICE_VENDOR,               sizeof d->vendor,    d->vendor,   NULL);
    clGetDeviceInfo(x, CL_DEVICE_VERSION,              sizeof d->clver,     d->clver,    NULL);
    clGetDeviceInfo(x, CL_DRIVER_VERSION,              sizeof d->driver,    d->driver,   NULL);
    clGetDeviceInfo(x, CL_DEVICE_GLOBAL_MEM_SIZE,      sizeof d->gmem,     &d->gmem,     NULL);
    clGetDeviceInfo(x, CL_DEVICE_LOCAL_MEM_SIZE,       sizeof d->lmem,     &d->lmem,     NULL);
    clGetDeviceInfo(x, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,sizeof d->cache,    &d->cache,    NULL);
    clGetDeviceInfo(x, CL_DEVICE_MAX_MEM_ALLOC_SIZE,   sizeof d->maxAlloc, &d->maxAlloc, NULL);
    clGetDeviceInfo(x, CL_DEVICE_MAX_COMPUTE_UNITS,    sizeof d->cus,      &d->cus,      NULL);
    clGetDeviceInfo(x, CL_DEVICE_MAX_CLOCK_FREQUENCY,  sizeof d->clockMHz, &d->clockMHz, NULL);
    clGetDeviceInfo(x, CL_DEVICE_MAX_WORK_GROUP_SIZE,  sizeof d->maxWG,    &d->maxWG,    NULL);
    d->isGpu = (d->type & CL_DEVICE_TYPE_GPU) != 0;
}

static void enumerateDevices(const char *want)
{
    cl_uint np = 0;
    if (clGetPlatformIDs(0, NULL, &np) != CL_SUCCESS || !np) return;
    cl_platform_id *plats = malloc((size_t)np * sizeof *plats);
    clGetPlatformIDs(np, plats, NULL);

    for (cl_uint i = 0; i < np && g_ndev < MAXDEV; i++) {
        cl_uint nd = 0;
        if (clGetDeviceIDs(plats[i], CL_DEVICE_TYPE_ALL, 0, NULL, &nd) != CL_SUCCESS || !nd) continue;
        cl_device_id *ds = malloc((size_t)nd * sizeof *ds);
        clGetDeviceIDs(plats[i], CL_DEVICE_TYPE_ALL, nd, ds, NULL);
        for (cl_uint j = 0; j < nd && g_ndev < MAXDEV; j++) {
            DevInfo t; memset(&t, 0, sizeof t);
            t.id = ds[j];
            collectDevice(&t);
            int isCpu = (t.type & CL_DEVICE_TYPE_CPU) != 0;
            if (!strcmp(want, "gpu") && !t.isGpu) continue;
            if (!strcmp(want, "cpu") && !isCpu)   continue;
            g_devs[g_ndev++] = t;
        }
        free(ds);
    }
    free(plats);

    if (getenv("MPA_FAKE_CPU_DEVICE") && g_ndev == 1 && g_ndev < MAXDEV) {
        g_devs[1] = g_devs[0];
        g_devs[1].isGpu = 0;
        g_devs[1].type = CL_DEVICE_TYPE_CPU;
        snprintf(g_devs[1].name, sizeof g_devs[1].name, "%s (posing as CPU)", g_devs[0].name);
        g_ndev = 2;
    }

    g_primary = 0;
    for (int i = 0; i < g_ndev; i++)
        if (g_devs[i].isGpu) { g_primary = i; break; }
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
    const char *wantDev = "all";

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--items")  && i+1 < argc) items  = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--reps")   && i+1 < argc) reps   = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--budget") && i+1 < argc) budget = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--variant")&& i+1 < argc) only   = argv[++i];
        else if (!strcmp(argv[i], "--devices")&& i+1 < argc) wantDev= argv[++i];
        else if (!strcmp(argv[i], "--verbose")) g_verbose = 1;
        else {
            fprintf(stderr,
                "usage: %s [--items N] [--reps N] [--budget SECONDS] "
                "[--variant w8|w16|w32|w32-opt|w32-o64] [--devices all|gpu|cpu] [--verbose]\n", argv[0]);
            return 2;
        }
    }
    if (items < 64) items = 64;
    if (reps  < 1)  reps  = 1;

    signal(SIGINT, onSigint);

    rngState = 88172645463325252ULL;

    enumerateDevices(wantDev);
    if (!g_ndev) { fprintf(stderr, "no OpenCL device matched --devices %s\n", wantDev); return 1; }

    detectCpu(&g_env);
    detectOs(&g_env);
#ifdef _OPENMP
    g_env.threads = omp_get_max_threads();
#else
    g_env.threads = 1;
#endif

    char safe[256];
    sanitize(g_devs[g_primary].name, safe, sizeof safe);
    snprintf(g_reportPath, sizeof g_reportPath, "%s_Report.md", safe);
    snprintf(g_csvPath,    sizeof g_csvPath,    "%s_Report.csv", safe);

    fprintf(stderr, "OpenCL devices under test (%d):\n", g_ndev);
    for (int i = 0; i < g_ndev; i++)
        fprintf(stderr, "  [%d] %-4s %s (%s), %u CU, %.2f GiB%s\n", i,
                g_devs[i].isGpu ? "GPU" : ((g_devs[i].type & CL_DEVICE_TYPE_CPU) ? "CPU" : "ACC"),
                g_devs[i].name, g_devs[i].vendor, g_devs[i].cus,
                (double)g_devs[i].gmem / 1073741824.0,
                i == g_primary ? "  <- names the report" : "");
    fprintf(stderr, "host CPU: %s, %ld cores, %.1f GB\n", g_env.cpu, g_env.cores, g_env.ramGB);
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

    for (int d = 0; d < g_ndev && !g_interrupted; d++) {
    cl_device_id dev = g_devs[d].id;
    cl_context ctx = clCreateContext(NULL, 1, &dev, NULL, NULL, &err);
    if (err != CL_SUCCESS) { fprintf(stderr, "context failed for %s\n", g_devs[d].name); continue; }
    cl_command_queue q = clCreateCommandQueue(ctx, dev, 0, &err);
    if (err != CL_SUCCESS) { clReleaseContext(ctx); continue; }
    fprintf(stderr, "\n=== device %d: %s ===\n", d, g_devs[d].name);

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
                fprintf(stderr, "  [%d] BUILD FAILED %s/%s: %s\n", d, var->name, mod->name, log);
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
                GpuCell *cell = &g_gpu[d][v][m][o];
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

                fprintf(stderr, "  [%d] %-8s %-18s %-26s n=%-7ld %.6fs %s\n",
                        d, var->name, modShort(mod), op->name, n, cell->kernel_s,
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

    clReleaseCommandQueue(q);
    clReleaseContext(ctx);
    }

    writeReport(items, reps, now_s() - t_start, !g_interrupted);

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

static int bestVariant(int d, int m, int o, double *secs)
{
    int best = -1; double b = 0;
    for (int v = 0; v < NVARIANTS; v++) {
        const GpuCell *c = &g_gpu[d][v][m][o];
        if (!c->built || c->mismatches || c->kernel_s <= 0) continue;
        if (best < 0 || c->kernel_s < b) { b = c->kernel_s; best = v; }
    }
    if (secs) *secs = b;
    return best;
}

static int bestOfClass(int wantGpu, int m, int o, double *secs, int *devOut, int *varOut)
{
    int found = 0; double b = 0; int bd = -1, bv = -1;
    for (int d = 0; d < g_ndev; d++) {
        if (g_devs[d].isGpu != wantGpu) continue;
        double s; int v = bestVariant(d, m, o, &s);
        if (v < 0) continue;
        if (!found || s < b) { b = s; bd = d; bv = v; found = 1; }
    }
    if (secs) *secs = b;
    if (devOut) *devOut = bd;
    if (varOut) *varOut = bv;
    return found;
}

static const char *devClass(const DevInfo *d)
{
    if (d->isGpu) return "GPU";
    if (d->type & CL_DEVICE_TYPE_CPU) return "CPU";
    return "ACC";
}

static void deviceTable(FILE *f, const DevInfo *d)
{
    fprintf(f, "| Property | Value |\n|---|---|\n");
    fprintf(f, "| Model | %s |\n", d->name);
    fprintf(f, "| Type | %s |\n", devClass(d));
    fprintf(f, "| Vendor | %s |\n", d->vendor);
    fprintf(f, "| Device memory | %.2f GiB |\n", (double)d->gmem / 1073741824.0);
    fprintf(f, "| Max single allocation | %.2f GiB |\n", (double)d->maxAlloc / 1073741824.0);
    fprintf(f, "| Local memory | %.0f KiB |\n", (double)d->lmem / 1024.0);
    fprintf(f, "| Global cache | %.0f KiB |\n", (double)d->cache / 1024.0);
    fprintf(f, "| Compute units | %u |\n", d->cus);
    fprintf(f, "| Max clock | %u MHz |\n", d->clockMHz);
    fprintf(f, "| Max work-group size | %zu |\n", d->maxWG);
    fprintf(f, "| OpenCL version | %s |\n", d->clver);
    fprintf(f, "| Driver | %s |\n\n", d->driver);
}

static void writeReport(int items, int reps, double elapsed, int complete)
{
    FILE *f = fopen(g_reportPath, "w");
    if (!f) { perror(g_reportPath); return; }

    fprintf(f, "# MPA-OpenCL benchmark report - %s\n\n", g_devs[g_primary].name);
    if (!complete)
        fprintf(f, "> **Partial report.** The run was interrupted or hit its time budget.\n"
                   "> Rows that never ran are marked `n/a`.\n\n");

    fprintf(f, "## 1. System under test\n\n");
    fprintf(f, "%d OpenCL device(s) exercised with the identical kernels and operands.\n\n", g_ndev);
    for (int d = 0; d < g_ndev; d++) {
        fprintf(f, "### Device %d - %s (%s)\n\n", d, g_devs[d].name, devClass(&g_devs[d]));
        deviceTable(f, &g_devs[d]);
    }

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
    fprintf(f, "- Every OpenCL device runs the same kernels on the same operands, so GPU and CPU-OpenCL columns are directly comparable.\n");
    fprintf(f, "- CPU library baselines (GMP, OpenSSL) run those same operands, with temporaries preallocated outside the timed region, so the figure is the arithmetic and not marshalling. The generator is reseeded per modulus and operation so every backend sees identical inputs.\n");
    fprintf(f, "- OpenSSL rows time the nearest BN primitive, which is not always semantically identical (its Montgomery routine expects Montgomery-domain inputs); they measure comparable work, not identical results. Correctness is judged against GMP only.\n");
    fprintf(f, "- Every device cell is checked word-for-word against GMP before it is timed. A cell that mismatches is reported and excluded from the speedup tables.\n");
    fprintf(f, "- Total wall time %.1f s.\n\n", elapsed);

    fprintf(f, "## 3. Correctness\n\n");
    fprintf(f, "| Device | Kernel | Configs run | Passed | Mismatched | Build/launch failed |\n|---|---|---|---|---|---|\n");
    long gtot = 0, gbad = 0;
    for (int d = 0; d < g_ndev; d++)
        for (int v = 0; v < NVARIANTS; v++) {
            long run = 0, pass = 0, bad = 0, fail = 0;
            for (int m = 0; m < NMODULI; m++)
                for (int o = 0; o < NOPS; o++) {
                    const GpuCell *c = &g_gpu[d][v][m][o];
                    if (!c->attempted) continue;
                    run++;
                    if (!c->built) fail++;
                    else if (c->mismatches) bad++;
                    else pass++;
                }
            if (!run) continue;
            gtot += run; gbad += bad + fail;
            fprintf(f, "| [%d] %s | `%s` (%s) | %ld | %ld | %ld | %ld |\n",
                    d, devClass(&g_devs[d]), VARIANTS[v].cl, VARIANTS[v].name, run, pass, bad, fail);
        }
    fprintf(f, "\n**%s** - %ld configurations, %ld problems.\n\n",
            gbad ? "FAILURES PRESENT" : "All configurations correct", gtot, gbad);

    fprintf(f, "## 4. Throughput per device\n\n");
    fprintf(f, "Operations per second, higher is better. Kernel-only timings.\n\n");
    for (int d = 0; d < g_ndev; d++) {
        fprintf(f, "### Device %d - %s (%s)\n\n", d, g_devs[d].name, devClass(&g_devs[d]));
        for (int m = 0; m < NMODULI; m++) {
            fprintf(f, "#### %s (%d-bit)\n\n| Operation | items |", MODULI[m].name, MODULI[m].bits);
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
                    const GpuCell *c = &g_gpu[d][v][m][o];
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
    }

    fprintf(f, "## 5. Head to head\n\n");
    fprintf(f, "Best OpenCL GPU result against best OpenCL CPU result and the CPU libraries.\n");
    fprintf(f, "Ratios above 1.00x mean the GPU is faster than that baseline.\n\n");
    for (int m = 0; m < NMODULI; m++) {
        fprintf(f, "### %s (%d-bit)\n\n", MODULI[m].name, MODULI[m].bits);
        fprintf(f, "| Operation | best GPU | GPU ops/s | best CPU-CL | CPU-CL ops/s | GMP 1T | GMP %dT | OpenSSL | CGBN |"
                   " GPU vs CPU-CL | GPU vs GMP %dT | GPU vs OpenSSL | GPU vs CGBN |\n",
                   g_env.threads, g_env.threads);
        fprintf(f, "|---|---|---|---|---|---|---|---|---|---|---|---|---|\n");
        for (int o = 0; o < NOPS; o++) {
            const CpuRow *cr = &g_cpu[m][o];
            if (!cr->items) continue;
            double gs = 0, cs2 = 0; int gd = -1, gv = -1, cd = -1, cv = -1;
            int haveG = bestOfClass(1, m, o, &gs,  &gd, &gv);
            int haveC = bestOfClass(0, m, o, &cs2, &cd, &cv);
            char gr[32], cr2[32], b1[32], b2[32], b3[32], b4[32];
            char r0[32], r1[32], r2[32], r3[32];
            long ci = 0; double cg = cgbnLookup(MODULI[m].name, OPS[o].name, &ci);
            double cgs = (cg > 0 && ci > 0) ? cg * (double)cr->items / (double)ci : -1;

            rate(gr,  sizeof gr,  haveG ? gs  : -1, cr->items);
            rate(cr2, sizeof cr2, haveC ? cs2 : -1, cr->items);
            rate(b1, sizeof b1, cr->gmp1, cr->items);
            rate(b2, sizeof b2, cr->gmpN, cr->items);
            rate(b3, sizeof b3, cr->ossl, cr->items);
            rate(b4, sizeof b4, cg, ci);
            ratio(r0, sizeof r0, haveC ? cs2 : -1, haveG ? gs : -1);
            ratio(r1, sizeof r1, cr->gmpN, haveG ? gs : -1);
            ratio(r2, sizeof r2, cr->ossl, haveG ? gs : -1);
            ratio(r3, sizeof r3, cgs,      haveG ? gs : -1);

            fprintf(f, "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n",
                    OPS[o].name,
                    haveG ? VARIANTS[gv].name : "none", gr,
                    haveC ? VARIANTS[cv].name : "none", cr2,
                    b1, b2, b3, b4, r0, r1, r2, r3);
        }
        fprintf(f, "\n");
    }

    int sec = 6;
    if (!g_ncgbn) {
        fprintf(f, "## %d. CGBN\n\n", sec++);
        fprintf(f, "Not measured on this run. CGBN is CUDA-only and is not built into this host.\n");
        fprintf(f, "To populate the CGBN columns, produce `cgbn_results.tsv` next to the binary with one\n");
        fprintf(f, "whitespace-separated row per measurement and re-run:\n\n");
        fprintf(f, "```\n# modulus_name  operation_name  items  seconds\nsecp256k1  MODMUL  20000  0.00123\n```\n\n");
        fprintf(f, "`modulus_name` and `operation_name` must match the spellings used in the tables above.\n\n");
    }

    fprintf(f, "## %d. Raw data\n\n", sec);
    fprintf(f, "Also written to `%s` for analysis.\n\n", g_csvPath);
    fprintf(f, "```csv\nkind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches\n");
    FILE *c = fopen(g_csvPath, "w");
    if (c) fprintf(c, "kind,device,device_type,kernel,modulus,bits,operation,items,seconds,ops_per_sec,mismatches\n");

    for (int m = 0; m < NMODULI; m++)
        for (int o = 0; o < NOPS; o++) {
            const CpuRow *cr = &g_cpu[m][o];
            if (!cr->items) continue;
            struct { const char *k; double s; } cpu[] = {
                { "gmp-1t", cr->gmp1 }, { "gmp-nt", cr->gmpN }, { "openssl-nt", cr->ossl } };
            for (int i = 0; i < 3; i++) {
                if (cpu[i].s <= 0) continue;
                const char *fmt = "library,%s,host-cpu,%s,%s,%d,%s,%ld,%.9f,%.3f,0\n";
                fprintf(f, fmt, g_env.cpu, cpu[i].k, MODULI[m].name, MODULI[m].bits,
                        OPS[o].name, cr->items, cpu[i].s, (double)cr->items / cpu[i].s);
                if (c) fprintf(c, fmt, g_env.cpu, cpu[i].k, MODULI[m].name, MODULI[m].bits,
                        OPS[o].name, cr->items, cpu[i].s, (double)cr->items / cpu[i].s);
            }
            long ci = 0; double cg = cgbnLookup(MODULI[m].name, OPS[o].name, &ci);
            if (cg > 0 && ci > 0) {
                const char *fmt = "library,%s,gpu,cgbn,%s,%d,%s,%ld,%.9f,%.3f,0\n";
                fprintf(f, fmt, g_devs[g_primary].name, MODULI[m].name, MODULI[m].bits,
                        OPS[o].name, ci, cg, (double)ci / cg);
                if (c) fprintf(c, fmt, g_devs[g_primary].name, MODULI[m].name, MODULI[m].bits,
                        OPS[o].name, ci, cg, (double)ci / cg);
            }
            for (int d = 0; d < g_ndev; d++)
                for (int v = 0; v < NVARIANTS; v++) {
                    const GpuCell *g = &g_gpu[d][v][m][o];
                    if (!g->attempted || !g->built) continue;
                    const char *fk = "opencl-kernel,%s,%s,%s,%s,%d,%s,%ld,%.9f,%.3f,%ld\n";
                    const char *fe = "opencl-e2e,%s,%s,%s,%s,%d,%s,%ld,%.9f,%.3f,%ld\n";
                    fprintf(f, fk, g_devs[d].name, devClass(&g_devs[d]), VARIANTS[v].name,
                            MODULI[m].name, MODULI[m].bits, OPS[o].name, g->items,
                            g->kernel_s, (double)g->items / g->kernel_s, g->mismatches);
                    fprintf(f, fe, g_devs[d].name, devClass(&g_devs[d]), VARIANTS[v].name,
                            MODULI[m].name, MODULI[m].bits, OPS[o].name, g->items,
                            g->e2e_s, (double)g->items / g->e2e_s, g->mismatches);
                    if (c) {
                        fprintf(c, fk, g_devs[d].name, devClass(&g_devs[d]), VARIANTS[v].name,
                                MODULI[m].name, MODULI[m].bits, OPS[o].name, g->items,
                                g->kernel_s, (double)g->items / g->kernel_s, g->mismatches);
                        fprintf(c, fe, g_devs[d].name, devClass(&g_devs[d]), VARIANTS[v].name,
                                MODULI[m].name, MODULI[m].bits, OPS[o].name, g->items,
                                g->e2e_s, (double)g->items / g->e2e_s, g->mismatches);
                    }
                }
        }
    fprintf(f, "```\n");
    if (c) fclose(c);
    fclose(f);
}
