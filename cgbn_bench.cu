// CGBN reference benchmark for MPA-OpenCL.
//
// Emits the cgbn_results.tsv that GPU_Host folds into its report.
//
//   ./GPU_Host --dump-moduli > moduli.tsv
//   ./cgbn_bench moduli.tsv [items] [reps] > cgbn_results.tsv
//
// Build:  make cgbn CGBN_PATH=/path/to/CGBN CUDA_ARCH=-arch=sm_90
// CGBN:   https://github.com/NVlabs/CGBN   (header only, needs nvcc and GMP)
//
// Operation names match GPU_Host's OPS table so the report can join on them.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>
#include <cuda_runtime.h>
#include <cgbn/cgbn.h>

#define CUDA_OK(call) do {                                                     \
    cudaError_t _e = (call);                                                   \
    if (_e != cudaSuccess) {                                                   \
        fprintf(stderr, "%s:%d %s\n", __FILE__, __LINE__, cudaGetErrorString(_e)); \
        exit(EXIT_FAILURE);                                                    \
    }                                                                          \
} while (0)

template<uint32_t tpi, uint32_t bits>
class mpa_params_t {
public:
    static const uint32_t TPB          = 0;
    static const uint32_t MAX_ROTATION = 4;
    static const uint32_t SHM_LIMIT    = 0;
    static const uint32_t CONSTANT_TIME= 0;
    static const uint32_t TPI          = tpi;
    static const uint32_t BITS         = bits;
    static const uint32_t WINDOW_BITS  = 5;
};

template<uint32_t bits>
struct instance_t {
    cgbn_mem_t<bits> a, b, m, r2, r;
};

enum OpId { OP_ADD, OP_SUB, OP_MUL, OP_CMP, OP_REM, OP_DIV,
            OP_MODMUL, OP_MONTMUL, OP_MODEXP,
            OP_ADDMOD, OP_SUBMOD, OP_MODMUL_R2, OP_COUNT };

static const char *OP_NAME[OP_COUNT] = {
    "ADD", "SUBTRACT", "MULTIPLYPRODUCTSCANNING", "COMPARE", "REDUCE", "DIVIDE",
    "MODMUL", "MONTGOMERYMULTIPLICATION", "MODEXP",
    "ADDMOD", "SUBTRACTMOD", "MODMUL_R2"
};


template<class params>
__global__ void mpa_cgbn_kernel(cgbn_error_report_t *report,
                                instance_t<params::BITS> *inst,
                                uint32_t count, int op, uint32_t np0)
{
    const int32_t idx = (blockIdx.x * blockDim.x + threadIdx.x) / params::TPI;
    if (idx >= (int32_t)count) return;

    typedef cgbn_context_t<params::TPI, params>   context_t;
    typedef cgbn_env_t<context_t, params::BITS>   env_t;
    typedef typename env_t::cgbn_t                bn_t;
    typedef typename env_t::cgbn_wide_t           bn_wide_t;

    context_t bn_context(cgbn_report_monitor, report, idx);
    env_t     bn_env(bn_context.template env<env_t>());

    bn_t a, b, m, r, t;
    bn_wide_t w;

    cgbn_load(bn_env, a, &(inst[idx].a));
    cgbn_load(bn_env, b, &(inst[idx].b));
    cgbn_load(bn_env, m, &(inst[idx].m));

    switch (op) {
    case OP_ADD: cgbn_add(bn_env, r, a, b); break;
    case OP_SUB: cgbn_sub(bn_env, r, a, b); break;
    case OP_MUL: cgbn_mul(bn_env, r, a, b); break;
    case OP_CMP: cgbn_set_ui32(bn_env, r, (uint32_t)(cgbn_compare(bn_env, a, b) + 1)); break;
    case OP_REM: cgbn_rem(bn_env, r, a, m); break;
    case OP_DIV:
        /* MPA's divModN defines q=r=0 for a zero divisor, and makeCase feeds
         * b=0 in two directed cases. CGBN raises an error instead, which used
         * to drop the whole DIVIDE row. */
        if (cgbn_equals_ui32(bn_env, b, 0)) cgbn_set_ui32(bn_env, r, 0);
        else                                cgbn_div_rem(bn_env, r, t, a, b);
        break;
    case OP_MODMUL:
        cgbn_mul_wide(bn_env, w, a, b);
        cgbn_rem_wide(bn_env, r, w, m);
        break;
    case OP_MONTMUL:
        cgbn_mont_mul(bn_env, r, a, b, m, np0);
        break;
    case OP_MODEXP:
        cgbn_modular_power(bn_env, r, a, b, m);
        break;
    case OP_ADDMOD: {
        const uint32_t c = cgbn_add(bn_env, r, a, b);
        if (c || cgbn_compare(bn_env, r, m) >= 0) cgbn_sub(bn_env, r, r, m);
        break;
    }
    case OP_SUBMOD:
        cgbn_sub(bn_env, r, a, b);
        if (cgbn_compare(bn_env, a, b) < 0) cgbn_add(bn_env, r, r, m);
        break;
    case OP_MODMUL_R2: {
        /* Mirrors op_modmul_r2: Mont(Mont(a,b), R^2) == a*b mod m, using the
         * host-supplied R^2 and np0 rather than an in-kernel domain transfer. */
        bn_t r2;
        cgbn_load(bn_env, r2, &(inst[idx].r2));
        cgbn_mont_mul(bn_env, t, a, b, m, np0);
        cgbn_mont_mul(bn_env, r, t, r2, m, np0);
        break;
    }
    default:
        cgbn_set_ui32(bn_env, r, 0);
        break;
    }

    cgbn_store(bn_env, &(inst[idx].r), r);
}

template<uint32_t bits>
static void toCgbn(cgbn_mem_t<bits> &out, const mpz_t v)
{
    memset(&out, 0, sizeof(out));
    size_t words = 0;
    if (mpz_sgn(v) > 0)
        mpz_export(out._limbs, &words, -1, sizeof(uint32_t), 0, 0, v);
}

static uint64_t rngState = 88172645463325252ULL;
static uint32_t rnd32(void)
{
    rngState ^= rngState << 13;
    rngState ^= rngState >> 7;
    rngState ^= rngState << 17;
    return (uint32_t)(rngState >> 32);
}

static const int OP_MODULAR[OP_COUNT] = { 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1 };

static void randBelow(mpz_t out, const mpz_t bound, int bits)
{
    mpz_set_ui(out, 0);
    for (int i = 0; i < (bits + 31) / 32; i++) {
        mpz_mul_2exp(out, out, 32);
        mpz_add_ui(out, out, rnd32());
    }
    mpz_mod(out, out, bound);
}

/* Mirrors makeCase() in mpa_ref.h so both benchmarks see identical inputs. */
static void makeCase(int j, int modular, const mpz_t p, int bits, mpz_t a, mpz_t b)
{
    mpz_t lim;
    mpz_init(lim);
    mpz_ui_pow_ui(lim, 2, (unsigned long)bits);

    switch (j) {
    case 0:  mpz_set_ui(a, 0); mpz_set_ui(b, 0); break;
    case 1:  mpz_set_ui(a, 1); mpz_set_ui(b, 0); break;
    case 2:  mpz_set_ui(a, 0); mpz_set_ui(b, 1); break;
    case 3:  mpz_set_ui(a, 1); mpz_set_ui(b, 1); break;
    case 4:  mpz_sub_ui(a, p, 1); mpz_set_ui(b, 1); break;
    case 5:  mpz_set_ui(a, 1); mpz_sub_ui(b, p, 1); break;
    case 6:  mpz_sub_ui(a, p, 1); mpz_sub_ui(b, p, 1); break;
    case 7:  mpz_sub_ui(a, p, 1); mpz_set(b, a); break;
    case 8:  mpz_sub_ui(a, lim, 1); mpz_set_ui(b, 1); break;
    case 9:  mpz_sub_ui(a, lim, 1); mpz_sub_ui(b, lim, 1); break;
    case 10: mpz_sub_ui(a, lim, 1); mpz_set(b, a); break;
    case 11: mpz_set_ui(a, 1); mpz_ui_pow_ui(b, 2, (unsigned long)(bits - 1)); break;
    case 12: mpz_sub_ui(a, p, 1); mpz_tdiv_q_ui(b, p, 2); break;
    case 13: mpz_tdiv_q_ui(a, p, 3); mpz_sub_ui(b, p, 2); break;
    case 18:
    case 19:
        mpz_set(a, p);
        if (j == 19) { mpz_mul_ui(a, a, 2); mpz_add_ui(a, a, 1); }
        mpz_set(b, p);
        break;
    case 16:
    case 17:
        mpz_ui_pow_ui(a, 2, (unsigned long)(bits / 2));
        mpz_sub_ui(a, a, 1);
        mpz_mul(a, a, a);
        if (j == 17) mpz_sub_ui(a, a, 1);
        mpz_set_ui(b, (j == 16) ? 3 : 2);
        break;
    default:
        if (modular) { randBelow(a, p, bits); randBelow(b, p, bits); }
        else         { randBelow(a, lim, bits); randBelow(b, lim, bits); }
        break;
    }

    if (modular) { mpz_mod(a, a, p); mpz_mod(b, b, p); }
    else         { mpz_mod(a, a, lim); mpz_mod(b, b, lim); }
    mpz_clear(lim);
}

/* np0 = -m^-1 mod 2^32, the same constant MPA passes as m_prime. */
static uint32_t mprime32(const mpz_t p)
{
    mpz_t base, inv, mp;
    uint32_t r;
    mpz_inits(base, inv, mp, NULL);
    mpz_ui_pow_ui(base, 2, 32);
    if (mpz_invert(inv, p, base) == 0) mpz_set_ui(inv, 1);
    mpz_sub(mp, base, inv);
    r = (uint32_t)mpz_get_ui(mp);
    mpz_clears(base, inv, mp, NULL);
    return r;
}

template<class params>
static void runModulus(const char *name, const mpz_t p, int items, int reps)
{
    const uint32_t BITS = params::BITS;
    typedef instance_t<BITS> inst_t;
    const uint32_t np0 = mprime32(p);

    mpz_t r2;
    mpz_init(r2);
    mpz_ui_pow_ui(r2, 2, (unsigned long)(2 * BITS));
    mpz_mod(r2, r2, p);

    inst_t *hInst = (inst_t *)calloc((size_t)items, sizeof(inst_t));
    if (!hInst) { fprintf(stderr, "out of memory\n"); exit(EXIT_FAILURE); }

    inst_t *dInst = NULL;
    CUDA_OK(cudaMalloc((void **)&dInst, (size_t)items * sizeof(inst_t)));

    cgbn_error_report_t *report = NULL;
    CUDA_OK(cgbn_error_report_alloc(&report));

    const int TPB = (params::TPB == 0) ? 128 : params::TPB;
    const int IPB = TPB / params::TPI;

    cudaEvent_t evStart, evStop;
    CUDA_OK(cudaEventCreate(&evStart));
    CUDA_OK(cudaEventCreate(&evStop));
    float best_ms = -1.0f;

    for (int op = 0; op < OP_COUNT; op++) {
        const int n = items;

        {   /* same operands GPU_Host used for this operator */
            mpz_t a, b;
            mpz_inits(a, b, NULL);
            rngState = 88172645463325252ULL;
            for (int i = 0; i < n; i++) {
                makeCase(i, OP_MODULAR[op], p, (int)BITS, a, b);
                toCgbn<BITS>(hInst[i].a, a);
                toCgbn<BITS>(hInst[i].b, b);
                toCgbn<BITS>(hInst[i].m, p);
                toCgbn<BITS>(hInst[i].r2, r2);
            }
            mpz_clears(a, b, NULL);
            if (cudaMemcpy(dInst, hInst, (size_t)n * sizeof(inst_t),
                           cudaMemcpyHostToDevice) != cudaSuccess) {
                fprintf(stderr, "# %s %s: operand upload failed\n", name, OP_NAME[op]);
                return;
            }
        }
        const int nblocks = (n + IPB - 1) / IPB;
        cudaError_t e = cudaSuccess;
        int bad = 0;

        for (int r = 0; r < reps + 2 && !bad; r++) {
            const int timed = (r >= 2);
            if (timed) cudaEventRecord(evStart);
            mpa_cgbn_kernel<params><<<nblocks, TPB>>>(report, dInst, (uint32_t)n, op, np0);
            if (timed) cudaEventRecord(evStop);
            e = cudaDeviceSynchronize();
            if (e != cudaSuccess) {
                fprintf(stderr, "# %s %s: %s\n", name, OP_NAME[op], cudaGetErrorString(e));
                bad = 1;
                break;
            }
            if (cgbn_error_report_check(report)) {
                fprintf(stderr, "# %s %s: CGBN reported an error\n", name, OP_NAME[op]);
                cgbn_error_report_reset(report);
                bad = 1;
                break;
            }
            if (timed) {
                float ms = 0.0f;
                if (cudaEventElapsedTime(&ms, evStart, evStop) == cudaSuccess)
                    if (best_ms < 0 || ms < best_ms) best_ms = ms;
            }
        }

        if (bad) {
            // A faulted kernel poisons the context; nothing after this can run.
            if (e != cudaSuccess && e != cudaErrorLaunchTimeout) {
                fprintf(stderr, "# %s: context unusable, abandoning this modulus\n", name);
                return;
            }
            best_ms = -1.0f;
            continue;
        }
        if (best_ms > 0) {
            printf("%s\t%s\t%d\t%.9f\n", name, OP_NAME[op], n, best_ms / 1000.0f);
            fflush(stdout);
        }
        best_ms = -1.0f;
    }

    mpz_clear(r2);
    CUDA_OK(cudaEventDestroy(evStart));
    CUDA_OK(cudaEventDestroy(evStop));
    CUDA_OK(cgbn_error_report_free(report));
    CUDA_OK(cudaFree(dInst));
    free(hInst);
}

static void dispatch(const char *name, int bits, const mpz_t p, int items, int reps)
{
    switch (bits) {
    case 256:  runModulus< mpa_params_t< 8,  256> >(name, p, items, reps); break;
    case 512:  runModulus< mpa_params_t< 8,  512> >(name, p, items, reps); break;
    case 1024: runModulus< mpa_params_t<16, 1024> >(name, p, items, reps); break;
    case 2048: runModulus< mpa_params_t<32, 2048> >(name, p, items, reps); break;
    default:
        fprintf(stderr, "# %s: %d-bit not instantiated, skipped\n", name, bits);
        break;
    }
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr,
            "usage: %s moduli.tsv [items] [reps] > cgbn_results.tsv\n"
            "       produce moduli.tsv with: ./GPU_Host --dump-moduli > moduli.tsv\n", argv[0]);
        return 2;
    }
    const int items = (argc > 2) ? atoi(argv[2]) : 20000;
    const int reps  = (argc > 3) ? atoi(argv[3]) : 5;

    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 1; }

    int dev = 0;
    cudaDeviceProp prop;
    CUDA_OK(cudaGetDevice(&dev));
    CUDA_OK(cudaGetDeviceProperties(&prop, dev));

    printf("# CGBN benchmark, device %s, items=%d reps=%d\n", prop.name, items, reps);
    printf("# modulus\toperation\titems\tseconds\n");

    char name[64], hex[1024];
    int  bits;
    while (fscanf(f, "%63s %d %1023s", name, &bits, hex) == 3) {
        mpz_t p;
        mpz_init(p);
        if (mpz_set_str(p, hex, 16) != 0) {
            fprintf(stderr, "# %s: bad hex, skipped\n", name);
            mpz_clear(p);
            continue;
        }
        dispatch(name, bits, p, items, reps);
        mpz_clear(p);
    }
    fclose(f);
    return 0;
}
