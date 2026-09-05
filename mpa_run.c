/*
 * mpa_run.c - minimal MPA-OpenCL host.
 *
 * Standard C99 plus an OpenCL 1.2 ICD. No GMP, no OpenSSL.
 *
 * Everything the kernel consumes is a plain word array. Only two values are
 * derived, and neither needs arbitrary precision:
 *
 *   m'      one machine word, by Hensel lifting  (mprime32)
 *   R^2 p   a modular doubling loop              (compute_r2)
 *
 * R^2 is read by MODMUL_R2 alone; MODMUL and MODEXP build the Montgomery
 * domain in-kernel via toMontN, and the remaining operators ignore it.
 *
 * Word order is big-endian: index 0 is most significant, index T-1 least.
 * That matches op_add in the kernel, which propagates carry from T-1 down.
 *
 * Usage:
 *   ./mpa_run <op> <p-hex> <a-hex> <b-hex> [<a-hex> <b-hex> ...]
 *   ./mpa_run <op> <p-hex> <a-hex> [<a-hex> ...]          (unary operators)
 *
 * The modulus sets the operand width: T = ceil(hex digits of p / 8) words.
 * Operators that take no modulus still need it as a width argument.
 *
 *   MPA_KERNEL=file.cl   override the kernel source (default: the opt kernel)
 *   MPA_BUILD="-D..."    extra build options appended after -DWORDLENGTH_T
 *
 * The buffer layout here is the kernel's default: item g owns the contiguous
 * words [g*T, g*T+T). Building with -DMPA_INTERLEAVED=1 switches the kernel to
 * a strided layout that this host does not mirror, so it is rejected rather
 * than answered wrongly.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#define ADD                      1
#define SUBTRACT                 2
#define ADDMOD                   3
#define SUBTRACTMOD              4
#define MULTIPLYOPERANDSCANNING   5
#define MULTIPLYPRODUCTSCANNING  6
#define MONTGOMERYMULTIPLICATION 7
#define COMPARE                  8
#define REDUCE                   9
#define MODMUL                  10
#define MODEXP                  11
#define EXPONENTIATION          12
#define DIVIDE                  13
#define ISQRT                   14
#define MODMUL_R2               15

/* name, opcode, 2T-word output, single operand, needs m', needs R^2 */
static const struct {
    const char *name;
    int op, wide, unary, mont, r2;
} OPS[] = {
    { "add",       ADD,                      0, 0, 0, 0 },
    { "sub",       SUBTRACT,                 0, 0, 0, 0 },
    { "addmod",    ADDMOD,                   0, 0, 0, 0 },
    { "submod",    SUBTRACTMOD,              0, 0, 0, 0 },
    { "mulop",     MULTIPLYOPERANDSCANNING,   1, 0, 0, 0 },
    { "mul",       MULTIPLYPRODUCTSCANNING,  1, 0, 0, 0 },
    { "montmul",   MONTGOMERYMULTIPLICATION, 0, 0, 1, 0 },
    { "compare",   COMPARE,                  0, 0, 0, 0 },
    { "reduce",    REDUCE,                   0, 1, 0, 0 },
    { "modmul",    MODMUL,                   0, 0, 1, 0 },
    { "modexp",    MODEXP,                   0, 0, 1, 0 },
    { "exp",       EXPONENTIATION,           0, 0, 0, 0 },
    { "div",       DIVIDE,                   1, 0, 0, 0 },
    { "isqrt",     ISQRT,                    0, 1, 0, 0 },
    { "modmul_r2", MODMUL_R2,                0, 0, 1, 1 },
};
#define NOPS ((int)(sizeof(OPS) / sizeof(OPS[0])))

#define CHECK(x) do {                                                       \
    cl_int _e = (x);                                                        \
    if (_e != CL_SUCCESS) {                                                 \
        fprintf(stderr, "%s:%d: OpenCL error %d\n", __FILE__, __LINE__,     \
                (int)_e);                                                   \
        exit(EXIT_FAILURE);                                                 \
    }                                                                       \
} while (0)

static void die(const char *msg)
{
    fprintf(stderr, "%s\n", msg);
    exit(EXIT_FAILURE);
}

/* ---------------------------------------------------------------- bignum --
 * Four small routines over big-endian word arrays. Between them they cover
 * everything the host has to know about multi-precision arithmetic.
 */

static int hexval(char c)
{
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

/* Strip an 0x prefix and leading zeros; returns the significant digit count. */
static const char *hex_body(const char *hex, size_t *len)
{
    size_t n = strlen(hex);
    if (n > 2 && hex[0] == '0' && (hex[1] == 'x' || hex[1] == 'X')) {
        hex += 2;
        n -= 2;
    }
    while (n > 1 && *hex == '0') { hex++; n--; }
    *len = n;
    return hex;
}

/* Right-aligned hex into T words. Returns -1 if it does not fit or is not hex. */
static int hex_to_words(const char *hex, uint32_t *w, int T)
{
    size_t n;
    const char *p = hex_body(hex, &n);
    if (n == 0 || n > (size_t)T * 8) return -1;

    memset(w, 0, (size_t)T * sizeof *w);
    int wi = T - 1;
    size_t pos = n;
    while (pos > 0) {
        size_t chunk = pos >= 8 ? 8 : pos;
        pos -= chunk;
        uint32_t v = 0;
        for (size_t k = 0; k < chunk; k++) {
            int d = hexval(p[pos + k]);
            if (d < 0) return -1;
            v = (v << 4) | (uint32_t)d;
        }
        w[wi--] = v;
    }
    return 0;
}

static void print_hex(const uint32_t *w, int n)
{
    int i = 0;
    while (i < n - 1 && w[i] == 0) i++;
    printf("%x", w[i]);
    for (i++; i < n; i++) printf("%08x", w[i]);
}

static int cmp_words(const uint32_t *a, const uint32_t *b, int n)
{
    for (int i = 0; i < n; i++)
        if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return 0;
}

/* a -= b, wrapping mod 2^(32n); the borrow out is deliberately discarded. */
static void sub_wrap(uint32_t *a, const uint32_t *b, int n)
{
    uint32_t borrow = 0;
    for (int i = n - 1; i >= 0; i--) {
        uint64_t d = (uint64_t)a[i] - b[i] - borrow;
        a[i] = (uint32_t)d;
        borrow = (d >> 63) & 1u;
    }
}

/* r = 2r mod p, given r < p. */
static void mod_double(uint32_t *r, const uint32_t *p, int n)
{
    uint32_t carry = 0;
    for (int i = n - 1; i >= 0; i--) {
        uint32_t hi = r[i] >> 31;
        r[i] = (r[i] << 1) | carry;
        carry = hi;
    }
    /* 2r < 2p, so a single conditional subtraction suffices. When the shift
     * overflowed, the true value is 2^(32n) + r and exceeds p by definition,
     * and the wrapped subtraction lands on the right T-word result. */
    if (carry || cmp_words(r, p, n) >= 0) sub_wrap(r, p, n);
}

/* r2 = 2^(2*32*T) mod p. This is the only place the host does real
 * multi-precision work, and only MODMUL_R2 needs the result. */
static void compute_r2(uint32_t *r2, const uint32_t *p, int T)
{
    memset(r2, 0, (size_t)T * sizeof *r2);
    r2[T - 1] = 1;
    if (cmp_words(r2, p, T) >= 0) sub_wrap(r2, p, T);
    for (int i = 0; i < 64 * T; i++) mod_double(r2, p, T);
}

/* m' = -p^-1 mod 2^32. An inverse mod 2^32 depends only on p mod 2^32, so
 * this is a single-word computation: four Newton steps from a seed correct
 * to 3 bits carry it to 3 -> 6 -> 12 -> 24 -> 48 bits. */
static uint32_t mprime32(uint32_t p0)
{
    uint32_t x = p0;                       /* p0 odd => correct mod 2^3 */
    for (int i = 0; i < 4; i++) x *= 2u - p0 * x;
    return (uint32_t)(-(int32_t)x);
}

/* ---------------------------------------------------------------- OpenCL -- */

static char *read_file(const char *path, size_t *len)
{
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(EXIT_FAILURE); }
    if (fseek(f, 0, SEEK_END) != 0) die("seek failed");
    long n = ftell(f);
    if (n < 0) die("ftell failed");
    rewind(f);
    char *buf = malloc((size_t)n + 1);
    if (!buf) die("out of memory");
    if (fread(buf, 1, (size_t)n, f) != (size_t)n) die("short read");
    buf[n] = '\0';
    fclose(f);
    *len = (size_t)n;
    return buf;
}

static cl_device_id pick_device(void)
{
    cl_uint nplat = 0;
    clGetPlatformIDs(0, NULL, &nplat);
    if (nplat == 0) die("no OpenCL platform found");

    cl_platform_id *plat = malloc(nplat * sizeof *plat);
    CHECK(clGetPlatformIDs(nplat, plat, NULL));

    cl_device_id dev = NULL;
    const cl_device_type order[] = { CL_DEVICE_TYPE_GPU, CL_DEVICE_TYPE_ALL };
    for (int t = 0; t < 2 && !dev; t++)
        for (cl_uint i = 0; i < nplat && !dev; i++)
            clGetDeviceIDs(plat[i], order[t], 1, &dev, NULL);

    free(plat);
    if (!dev) die("no OpenCL device found");
    return dev;
}

int main(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr,
                "usage: %s <op> <p-hex> <a-hex> <b-hex> [<a-hex> <b-hex> ...]\n"
                "       %s <op> <p-hex> <a-hex> [<a-hex> ...]   (reduce, isqrt)\n"
                "ops:", argv[0], argv[0]);
        for (int i = 0; i < NOPS; i++) fprintf(stderr, " %s", OPS[i].name);
        fprintf(stderr, "\n");
        return EXIT_FAILURE;
    }

    int oi = -1;
    for (int i = 0; i < NOPS; i++)
        if (strcmp(argv[1], OPS[i].name) == 0) { oi = i; break; }
    if (oi < 0) die("unknown operator");

    /* The modulus fixes the width. */
    size_t plen;
    hex_body(argv[2], &plen);
    const int T = (int)((plen + 7) / 8);
    if (T < 1) die("modulus is empty");

    const int operands = argc - 3;
    const int per_item = OPS[oi].unary ? 1 : 2;
    if (operands % per_item != 0)
        die(OPS[oi].unary ? "expected one operand per item"
                          : "operands must come in pairs");
    const size_t items = (size_t)operands / per_item;

    const int outWords = OPS[oi].wide ? 2 * T : T;

    uint32_t *hA = calloc(items * (size_t)T, sizeof *hA);
    uint32_t *hB = calloc(items * (size_t)T, sizeof *hB);
    uint32_t *hC = calloc(items * (size_t)outWords, sizeof *hC);
    uint32_t *hP = calloc((size_t)T * 2, sizeof *hP);   /* p, then R^2 */
    if (!hA || !hB || !hC || !hP) die("out of memory");

    if (hex_to_words(argv[2], hP, T) != 0) die("bad modulus");

    for (size_t i = 0; i < items; i++) {
        const char *a = argv[3 + i * (size_t)per_item];
        if (hex_to_words(a, &hA[i * (size_t)T], T) != 0)
            die("operand does not fit in the modulus width");
        if (per_item == 2) {
            const char *b = argv[4 + i * (size_t)per_item];
            if (hex_to_words(b, &hB[i * (size_t)T], T) != 0)
                die("operand does not fit in the modulus width");
        }
    }

    uint32_t m_prime = 0;
    if (OPS[oi].mont) {
        if ((hP[T - 1] & 1u) == 0) die("Montgomery operators need an odd modulus");
        m_prime = mprime32(hP[T - 1]);
    }
    if (OPS[oi].r2) compute_r2(&hP[T], hP, T);

    cl_int opBuf[4];
    opBuf[0] = OPS[oi].op;
    opBuf[1] = 32;                    /* word size  */
    opBuf[2] = T * 32;                /* bit length */
    opBuf[3] = (cl_int)m_prime;

    cl_device_id dev = pick_device();
    cl_int err;
    cl_context ctx = clCreateContext(NULL, 1, &dev, NULL, NULL, &err); CHECK(err);
    cl_command_queue q = clCreateCommandQueue(ctx, dev, 0, &err); CHECK(err);

    const char *ksrc = getenv("MPA_KERNEL");
    if (!ksrc) ksrc = "mpaKernel_32bits_opt.cl";
    size_t srcLen;
    char *src = read_file(ksrc, &srcLen);

    const char *extra = getenv("MPA_BUILD");
    if (extra && strstr(extra, "MPA_INTERLEAVED"))
        die("MPA_INTERLEAVED changes the buffer layout; this host assumes the "
            "contiguous default");
    char opts[512];
    snprintf(opts, sizeof opts, "-I. -DWORDLENGTH_T=%d %s", T, extra ? extra : "");

    cl_program prog = clCreateProgramWithSource(ctx, 1, (const char **)&src,
                                                &srcLen, &err); CHECK(err);
    if (clBuildProgram(prog, 1, &dev, opts, NULL, NULL) != CL_SUCCESS) {
        size_t n = 0;
        clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &n);
        char *log = malloc(n + 1);
        clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, n, log, NULL);
        log[n] = '\0';
        fprintf(stderr, "build failed (%s):\n%s\n", opts, log);
        return EXIT_FAILURE;
    }
    free(src);

    cl_kernel kern = clCreateKernel(prog, "mpaKernel", &err); CHECK(err);

    const size_t inBytes  = items * (size_t)T * sizeof(uint32_t);
    const size_t outBytes = items * (size_t)outWords * sizeof(uint32_t);
    const size_t pBytes   = (size_t)T * 2 * sizeof(uint32_t);

    cl_mem dA = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  inBytes,      NULL, &err); CHECK(err);
    cl_mem dB = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  inBytes,      NULL, &err); CHECK(err);
    cl_mem dC = clCreateBuffer(ctx, CL_MEM_READ_WRITE, outBytes,     NULL, &err); CHECK(err);
    cl_mem dO = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  sizeof opBuf, NULL, &err); CHECK(err);
    cl_mem dP = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  pBytes,       NULL, &err); CHECK(err);

    CHECK(clSetKernelArg(kern, 0, sizeof(cl_mem), &dA));
    CHECK(clSetKernelArg(kern, 1, sizeof(cl_mem), &dB));
    CHECK(clSetKernelArg(kern, 2, sizeof(cl_mem), &dC));
    CHECK(clSetKernelArg(kern, 3, sizeof(cl_mem), &dO));
    CHECK(clSetKernelArg(kern, 4, sizeof(cl_mem), &dP));

    CHECK(clEnqueueWriteBuffer(q, dO, CL_TRUE, 0, sizeof opBuf, opBuf, 0, NULL, NULL));
    CHECK(clEnqueueWriteBuffer(q, dP, CL_TRUE, 0, pBytes,  hP, 0, NULL, NULL));
    CHECK(clEnqueueWriteBuffer(q, dA, CL_TRUE, 0, inBytes, hA, 0, NULL, NULL));
    CHECK(clEnqueueWriteBuffer(q, dB, CL_TRUE, 0, inBytes, hB, 0, NULL, NULL));

    CHECK(clEnqueueNDRangeKernel(q, kern, 1, NULL, &items, NULL, 0, NULL, NULL));
    CHECK(clFinish(q));
    CHECK(clEnqueueReadBuffer(q, dC, CL_TRUE, 0, outBytes, hC, 0, NULL, NULL));

    for (size_t i = 0; i < items; i++) {
        print_hex(&hC[i * (size_t)outWords], outWords);
        printf("\n");
    }

    clReleaseMemObject(dA); clReleaseMemObject(dB); clReleaseMemObject(dC);
    clReleaseMemObject(dO); clReleaseMemObject(dP);
    clReleaseKernel(kern);
    clReleaseProgram(prog);
    clReleaseCommandQueue(q);
    clReleaseContext(ctx);
    free(hA); free(hB); free(hC); free(hP);
    return EXIT_SUCCESS;
}
