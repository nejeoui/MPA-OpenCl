#ifndef MPA_REF_H
#define MPA_REF_H

#define CL_TARGET_OPENCL_VERSION 120
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#define ADD                      1
#define SUBTRACT                 2
#define ADDMOD                   3
#define SUBTRACTMOD              4
#define MULTIPLYOPRANDSCANNING   5
#define MULTIPLYPRODUCTSCANNING  6
#define MONTGOMERYMULTIPLICATION 7
#define COMPARE 8
#define REDUCE 9
#define MODMUL 10
#define MODEXP 11
#define EXPONENTIATION 12
#define DIVIDE 13
#define ISQRT 14
#define MODMUL_R2 15

#define MAX_SOURCE_SIZE (0x100000)

#define GREEN "\x1b[32m"
#define RED   "\x1b[31m"
#define DIM   "\x1b[2m"
#define OFF   "\x1b[0m"

static int g_verbose = 0;

static const char *clErr(cl_int e)
{
    switch (e) {
    case CL_SUCCESS:                    return "CL_SUCCESS";
    case CL_DEVICE_NOT_FOUND:           return "CL_DEVICE_NOT_FOUND";
    case CL_BUILD_PROGRAM_FAILURE:      return "CL_BUILD_PROGRAM_FAILURE";
    case CL_INVALID_VALUE:              return "CL_INVALID_VALUE";
    case CL_INVALID_PLATFORM:           return "CL_INVALID_PLATFORM";
    case CL_INVALID_DEVICE:             return "CL_INVALID_DEVICE";
    case CL_INVALID_KERNEL_NAME:        return "CL_INVALID_KERNEL_NAME";
    case CL_INVALID_KERNEL_ARGS:        return "CL_INVALID_KERNEL_ARGS";
    case CL_INVALID_WORK_GROUP_SIZE:    return "CL_INVALID_WORK_GROUP_SIZE";
    case CL_MEM_OBJECT_ALLOCATION_FAILURE: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
    case CL_OUT_OF_RESOURCES:           return "CL_OUT_OF_RESOURCES";
    case CL_OUT_OF_HOST_MEMORY:         return "CL_OUT_OF_HOST_MEMORY";
    default:                            return "CL error";
    }
}

#define CHECK(expr) do {                                                     \
    cl_int _e = (expr);                                                      \
    if (_e != CL_SUCCESS) {                                                  \
        fprintf(stderr, "%s:%d: %s -> %s (%d)\n",                            \
                __FILE__, __LINE__, #expr, clErr(_e), (int)_e);              \
        exit(2);                                                             \
    }                                                                        \
} while (0)

static const char *devTypeName(cl_device_type t)
{
    if (t & CL_DEVICE_TYPE_GPU)         return "GPU";
    if (t & CL_DEVICE_TYPE_ACCELERATOR) return "ACCELERATOR";
    if (t & CL_DEVICE_TYPE_CPU)         return "CPU";
    return "OTHER";
}

static void listDevices(void)
{
    cl_platform_id plats[16];
    cl_uint nplat = 0;
    if (clGetPlatformIDs(16, plats, &nplat) != CL_SUCCESS) return;
    fprintf(stderr, "available OpenCL devices:\n");
    for (cl_uint i = 0; i < nplat; i++) {
        cl_device_id devs[16];
        cl_uint nd = 0;
        char pname[256] = {0};
        clGetPlatformInfo(plats[i], CL_PLATFORM_NAME, sizeof(pname), pname, NULL);
        if (clGetDeviceIDs(plats[i], CL_DEVICE_TYPE_ALL, 16, devs, &nd) != CL_SUCCESS) continue;
        for (cl_uint j = 0; j < nd; j++) {
            char dname[256] = {0};
            cl_device_type dt = 0;
            clGetDeviceInfo(devs[j], CL_DEVICE_NAME, sizeof(dname), dname, NULL);
            clGetDeviceInfo(devs[j], CL_DEVICE_TYPE, sizeof(dt), &dt, NULL);
            fprintf(stderr, "  [%u.%u] %-11s %s (%s)\n", i, j, devTypeName(dt), dname, pname);
        }
    }
}

static void pickDevice(cl_platform_id *outPlat, cl_device_id *outDev)
{
    cl_platform_id plats[16];
    cl_uint nplat = 0;
    const cl_device_type order[3] = { CL_DEVICE_TYPE_GPU,
                                      CL_DEVICE_TYPE_ACCELERATOR,
                                      CL_DEVICE_TYPE_CPU };
    const char *want = getenv("MPA_DEVICE_TYPE");
    const char *idxs = getenv("MPA_DEVICE_INDEX");
    long wantIdx = idxs ? strtol(idxs, NULL, 10) : 0;
    int t0 = 0, t1 = 3, strict = 0;

    if (getenv("MPA_LIST_DEVICES")) { listDevices(); exit(0); }

    CHECK(clGetPlatformIDs(16, plats, &nplat));
    if (nplat == 0) { fprintf(stderr, "no OpenCL platform found\n"); exit(2); }

    if (want) {
        strict = 1;
        if      (!strcmp(want, "gpu"))         { t0 = 0; t1 = 1; }
        else if (!strcmp(want, "accelerator")) { t0 = 1; t1 = 2; }
        else if (!strcmp(want, "cpu"))         { t0 = 2; t1 = 3; }
        else if (!strcmp(want, "any"))         { t0 = 0; t1 = 3; strict = 0; }
        else {
            fprintf(stderr, "MPA_DEVICE_TYPE must be gpu, accelerator, cpu or any\n");
            exit(2);
        }
    }

    long seen = 0;
    for (int t = t0; t < t1; t++) {
        for (cl_uint i = 0; i < nplat; i++) {
            cl_device_id devs[16];
            cl_uint nd = 0;
            if (clGetDeviceIDs(plats[i], order[t], 16, devs, &nd) != CL_SUCCESS) continue;
            for (cl_uint j = 0; j < nd; j++) {
                if (seen++ != wantIdx) continue;
                *outPlat = plats[i];
                *outDev  = devs[j];
                return;
            }
        }
    }

    if (strict) {
        fprintf(stderr,
                "no %s device found (MPA_DEVICE_TYPE=%s, MPA_DEVICE_INDEX=%ld).\n"
                "Refusing to silently fall back to another device type.\n"
                "Run with MPA_LIST_DEVICES=1 to see what this machine exposes.\n",
                want, want, wantIdx);
    } else {
        fprintf(stderr, "no OpenCL device found\n");
    }
    exit(2);
}

static char *readFile(const char *path, size_t *len)
{
    FILE *fp = fopen(path, "rb");
    char *buf;
    if (!fp) { fprintf(stderr, "cannot open %s\n", path); exit(2); }
    buf = malloc(MAX_SOURCE_SIZE);
    if (!buf) { fprintf(stderr, "oom\n"); exit(2); }
    *len = fread(buf, 1, MAX_SOURCE_SIZE, fp);
    fclose(fp);
    return buf;
}

static void storeWord(void *buf, size_t idx, int wbits, uint32_t v)
{
    if      (wbits == 8)  ((uint8_t  *)buf)[idx] = (uint8_t)v;
    else if (wbits == 16) ((uint16_t *)buf)[idx] = (uint16_t)v;
    else                  ((uint32_t *)buf)[idx] = v;
}

static uint32_t loadWord(const void *buf, size_t idx, int wbits)
{
    if      (wbits == 8)  return ((const uint8_t  *)buf)[idx];
    else if (wbits == 16) return ((const uint16_t *)buf)[idx];
    else                  return ((const uint32_t *)buf)[idx];
}

static void mpzToWords(const mpz_t z, void *buf, size_t off, int nwords, int wbits)
{
    mpz_t t;
    uint32_t mask = (wbits == 32) ? 0xFFFFFFFFu : ((1u << wbits) - 1u);
    mpz_init_set(t, z);
    for (int i = nwords - 1; i >= 0; i--) {
        storeWord(buf, off + (size_t)i, wbits, (uint32_t)(mpz_get_ui(t) & mask));
        mpz_tdiv_q_2exp(t, t, (mp_bitcnt_t)wbits);
    }
    mpz_clear(t);
}

static void wordsToMpz(mpz_t z, const void *buf, size_t off, int nwords, int wbits)
{
    mpz_set_ui(z, 0);
    for (int i = 0; i < nwords; i++) {
        mpz_mul_2exp(z, z, (mp_bitcnt_t)wbits);
        mpz_add_ui(z, z, (unsigned long)loadWord(buf, off + (size_t)i, wbits));
    }
}

typedef struct { const char *name; int bits; const char *hex;
                 const char *facA, *facB; } Modulus;

static const Modulus MODULI[] = {
    { "secp256k1", 256,
      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
      NULL, NULL },
    { "rsa256(composite)", 256,
      "A9FC2F3C7C6E33B082ACB3C97EF3FEEEFDC744963673BDE28707F902982B585B",
      "B210C163FFCA64C508639DB7109A3467",
      "F46211F7526FD128BD09029B809363ED" },
    { "brainpoolP512r1", 512,
      "AADD9DB8DBE9C48B3FD4E6AE33C9FC07CB308DB3B3C9D20ED6639CCA70330871"
      "7D4D9B009BC66842AECDA12AE6A380E62881FF2F2D82C68528AA6056583A48F3",
      NULL, NULL },
    { "p1024", 1024,
      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF97",
      NULL, NULL },
    { "p2048", 2048,
      "E53DF5FC3F650D066875837012A4E7BEA863C65CB592D9C36942CF69CBC6DD4F"
      "D804E19CCF2696C9BEBCF18742FA5FB091CBDE1782E8291009464913ECE19745"
      "7800EA6E43B0E2A64615D182B6DE150479C58D1C7C702D47EA3031B379CA13A2"
      "048C964E1D1E8D4CD3815D0895BF31E53271D4607E16461B77FB26100915D679"
      "9060203EDEBFEA9495A5A8E7CED68FC9DB2D47CE7992461BA78174608AD0BBE3"
      "F5E63EC6C960564430CBD2E6E587D08EE12F94B5B99DFFB12C6727A25E800DAC"
      "6CD8DE77A5BBC93B36E444B070888CB5ADD991870466968A6E9A23C2EE0A1D67"
      "1C9B601081A44AA6A58D4DC76686EF15FCE1C9AEB4033395A9B24BE1AA1929BB",
      NULL, NULL },
};
#define NMODULI ((int)(sizeof(MODULI)/sizeof(MODULI[0])))

typedef struct { const char *cl; const char *name; int wbits; const char *flags; int ext; } Variant;

#define MPA_OPT_FLAGS "-DMPA_MULHI=1 -DMPA_REGACC=1 -DMPA_FUSED_CIOS=1 -DMPA_UNROLL=1"

static const Variant VARIANTS[] = {
    { "mpaKernels_8bits.cl",     "w8",      8,  "", 0 },
    { "mpaKernel_16bits.cl",     "w16",     16, "", 0 },
    { "mpaKernel_32bits.cl",     "w32",     32, "", 0 },
    { "mpaKernel_32bits_opt.cl", "w32-opt", 32, MPA_OPT_FLAGS, 1 },
    { "mpaKernel_32bits_opt.cl", "w32-o64", 32, "-DMPA_REGACC=1 -DMPA_FUSED_CIOS=1 -DMPA_UNROLL=1", 1 },
};
#define NVARIANTS ((int)(sizeof(VARIANTS)/sizeof(VARIANTS[0])))

typedef struct { int op; const char *name; int wide; int modular; int ext; int cost; } Op;

static const Op OPS[] = {
    { ADD,                      "ADD",                     0, 0, 0,  1 },
    { SUBTRACT,                 "SUBTRACT",                0, 0, 0,  1 },
    { ADDMOD,                   "ADDMOD",                  0, 1, 0,  1 },
    { SUBTRACTMOD,              "SUBTRACTMOD",             0, 1, 0,  1 },
    { MULTIPLYOPRANDSCANNING,   "MULTIPLYOPERANDSCANNING", 1, 0, 0,  1 },
    { MULTIPLYPRODUCTSCANNING,  "MULTIPLYPRODUCTSCANNING", 1, 0, 0,  1 },
    { MONTGOMERYMULTIPLICATION, "MONTGOMERYMULTIPLICATION",0, 1, 0,  1 },
    { COMPARE,                  "COMPARE",                 0, 0, 1,  1 },
    { REDUCE,                   "REDUCE",                  0, 0, 1,  8 },
    { MODMUL,                   "MODMUL",                  0, 1, 1, 16 },
    { MODEXP,                   "MODEXP",                  0, 1, 1, 64 },
    { EXPONENTIATION,           "EXPONENTIATION",          0, 0, 1, 64 },
    { DIVIDE,                   "DIVIDE",                  1, 0, 1,  8 },
    { ISQRT,                    "ISQRT",                   0, 0, 1, 32 },
    { MODMUL_R2,                "MODMUL_R2",               0, 1, 1,  1 },
};
#define NOPS ((int)(sizeof(OPS)/sizeof(OPS[0])))

#define EDGE_CASES 20

static uint64_t rngState;
static uint32_t rnd32(void)
{
    uint64_t z = (rngState += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return (uint32_t)((z ^ (z >> 31)) >> 16);
}

static void randBelow(mpz_t out, const mpz_t bound, int bits)
{
    mpz_t r;
    mpz_init_set_ui(r, 0);
    for (int i = 0; i < bits; i += 32) {
        mpz_mul_2exp(r, r, 32);
        mpz_add_ui(r, r, rnd32());
    }
    mpz_mod(out, r, bound);
    mpz_clear(r);
}

static void makeCase(int j, int op, int modular, const mpz_t p, int bits,
                     const Modulus *mod, mpz_t a, mpz_t b)
{
    mpz_t lim, one;
    mpz_init(lim);
    mpz_init_set_ui(one, 1);
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
    case 14:
    case 15:
        if (mod->facA && mod->facB) {
            mpz_set_str(a, j == 14 ? mod->facA : mod->facB, 16);
            mpz_set_str(b, j == 14 ? mod->facB : mod->facA, 16);
        } else {
            randBelow(a, modular ? p : p, bits);
            randBelow(b, modular ? p : p, bits);
        }
        break;
    default:
        if (modular) { randBelow(a, p, bits); randBelow(b, p, bits); }
        else         { randBelow(a, lim, bits); randBelow(b, lim, bits); }
        break;
    }

    if (modular) { mpz_mod(a, a, p); mpz_mod(b, b, p); }
    else         { mpz_mod(a, a, lim); mpz_mod(b, b, lim); }

    if (op == SUBTRACT && mpz_cmp(a, b) < 0) mpz_swap(a, b);

    mpz_clear(lim);
    mpz_clear(one);
}

static void computeR2(mpz_t r2, const mpz_t p, int bits)
{
    mpz_ui_pow_ui(r2, 2, (unsigned long)(2 * bits));
    mpz_mod(r2, r2, p);
}

static void reference(int op, const mpz_t a, const mpz_t b, const mpz_t p,
                      int bits, mpz_t out)
{
    mpz_t lim, R, Rinv;
    mpz_init(lim);
    mpz_ui_pow_ui(lim, 2, (unsigned long)bits);

    switch (op) {
    case ADD:
        mpz_add(out, a, b);
        mpz_mod(out, out, lim);
        break;
    case SUBTRACT:
        mpz_sub(out, a, b);
        mpz_mod(out, out, lim);
        break;
    case ADDMOD:
        mpz_add(out, a, b);
        mpz_mod(out, out, p);
        break;
    case SUBTRACTMOD:
        mpz_sub(out, a, b);
        mpz_mod(out, out, p);
        break;
    case MULTIPLYOPRANDSCANNING:
    case MULTIPLYPRODUCTSCANNING:
        mpz_mul(out, a, b);
        break;
    case MONTGOMERYMULTIPLICATION:
        mpz_init_set(R, lim);
        mpz_init(Rinv);
        mpz_invert(Rinv, R, p);
        mpz_mul(out, a, b);
        mpz_mul(out, out, Rinv);
        mpz_mod(out, out, p);
        mpz_clear(R);
        mpz_clear(Rinv);
        break;
    case COMPARE: {
        int c = mpz_cmp(a, b);
        c = (c > 0) - (c < 0);
        mpz_set_si(out, (long)c);
        mpz_mod(out, out, lim);
        break;
    }
    case REDUCE:
        mpz_mod(out, a, p);
        break;
    case MODMUL:
    case MODMUL_R2:
        mpz_mul(out, a, b);
        mpz_mod(out, out, p);
        break;
    case MODEXP:
        mpz_powm(out, a, b, p);
        break;
    case EXPONENTIATION:
        mpz_powm(out, a, b, lim);
        break;
    case DIVIDE: {
        mpz_t q, r;
        mpz_inits(q, r, NULL);
        if (mpz_sgn(b) != 0) mpz_tdiv_qr(q, r, a, b);
        mpz_mul_2exp(out, q, (mp_bitcnt_t)bits);
        mpz_add(out, out, r);
        mpz_clears(q, r, NULL);
        break;
    }
    case ISQRT:
        mpz_sqrt(out, a);
        break;
    default:
        mpz_set_ui(out, 0);
        break;
    }
    mpz_clear(lim);
}

#endif
