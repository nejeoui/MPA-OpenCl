#ifndef WORDLENGTH_T
#define WORDLENGTH_T 8
#endif

#ifndef MPA_MULHI
#define MPA_MULHI 0
#endif
#ifndef MPA_REGACC
#define MPA_REGACC 0
#endif
#ifndef MPA_FUSED_CIOS
#define MPA_FUSED_CIOS 0
#endif
#ifndef MPA_UNROLL
#define MPA_UNROLL 0
#endif
#ifndef MPA_INTERLEAVED
#define MPA_INTERLEAVED 0
#endif

#define ADD 1
#define SUBTRACT 2
#define ADDMOD 3
#define SUBTRACTMOD 4
#define MULTIPLYOPERANDSCANNING 5
#define MULTIPLYPRODUCTSCANNING 6
#define MONTGOMERYMULTIPLICATION 7
#define COMPARE 8
#define REDUCE 9
#define MODMUL 10
#define MODEXP 11
#define EXPONENTIATION 12
#define DIVIDE 13
#define ISQRT 14
#define MODMUL_R2 15

#define T WORDLENGTH_T

#if MPA_UNROLL
#define UNROLL __attribute__((opencl_unroll_hint))
#else
#define UNROLL
#endif

#if MPA_INTERLEAVED
  #define IDX_N(g,i,n)  ((size_t)(i) * get_global_size(0) + (size_t)(g))
#else
  #define IDX_N(g,i,n)  ((size_t)(g) * (size_t)(n) + (size_t)(i))
#endif
#define IDX(g,i)    IDX_N(g,i,T)
#define IDX2(g,i)   IDX_N(g,i,2*T)

inline uint addc(uint a, uint b, uint cin, uint *cout)
{
#if MPA_MULHI
    uint s = a + b;
    uint c = (s < a) ? 1u : 0u;
    s += cin;
    c += (s < cin) ? 1u : 0u;
    *cout = c;
    return s;
#else
    ulong s = (ulong)a + (ulong)b + (ulong)cin;
    *cout = (uint)(s >> 32);
    return (uint)s;
#endif
}

inline uint subb(uint a, uint b, uint bin, uint *bout)
{
#if MPA_MULHI
    uint d = a - b;
    uint br = (a < b) ? 1u : 0u;
    uint d2 = d - bin;
    br += (d < bin) ? 1u : 0u;
    *bout = br;
    return d2;
#else
    ulong d = (ulong)a - (ulong)b - (ulong)bin;
    *bout = (uint)((d >> 32) & 1u);
    return (uint)d;
#endif
}

inline void mac(uint a, uint b, uint c, uint d, uint *hi, uint *lo)
{
#if MPA_MULHI
    uint l = a * b;
    uint h = mul_hi(a, b);
    l += c; if (l < c) h++;
    l += d; if (l < d) h++;
    *hi = h; *lo = l;
#else
    ulong p = (ulong)a * (ulong)b + (ulong)c + (ulong)d;
    *hi = (uint)(p >> 32);
    *lo = (uint)p;
#endif
}

inline int cmpT(const uint *a, const uint *b)
{
    UNROLL
    for (int i = 0; i < T; i++) {
        if (a[i] > b[i]) return  1;
        if (a[i] < b[i]) return -1;
    }
    return 0;
}

inline void subInPlace(uint *r, const uint *p)
{
    uint borrow = 0;
    UNROLL
    for (int i = T - 1; i >= 0; i--) r[i] = subb(r[i], p[i], borrow, &borrow);
}

inline void addInPlace(uint *r, const uint *p)
{
    uint carry = 0;
    UNROLL
    for (int i = T - 1; i >= 0; i--) r[i] = addc(r[i], p[i], carry, &carry);
}

inline void loadT(uint *dst, __global const uint *src, size_t g)
{
    UNROLL
    for (int i = 0; i < T; i++) dst[i] = src[IDX(g, i)];
}

inline void storeT(__global uint *dst, const uint *src, size_t g)
{
    UNROLL
    for (int i = 0; i < T; i++) dst[IDX(g, i)] = src[i];
}

inline void op_add(__global const uint *x, __global const uint *y,
                   __global uint *out, size_t g)
{
    uint carry = 0;
    UNROLL
    for (int i = T - 1; i >= 0; i--) {
        const size_t k = IDX(g, i);
        out[k] = addc(x[k], y[k], carry, &carry);
    }
}

inline void op_sub(__global const uint *x, __global const uint *y,
                   __global uint *out, size_t g)
{
    uint borrow = 0;
    UNROLL
    for (int i = T - 1; i >= 0; i--) {
        const size_t k = IDX(g, i);
        out[k] = subb(x[k], y[k], borrow, &borrow);
    }
}

inline void op_addmod(__global const uint *x, __global const uint *y,
                      __global uint *out, size_t g, const uint *p)
{
#if MPA_REGACC
    uint r[T];
    uint carry = 0;
    UNROLL
    for (int i = T - 1; i >= 0; i--) {
        const size_t k = IDX(g, i);
        r[i] = addc(x[k], y[k], carry, &carry);
    }
    if (carry || cmpT(r, p) >= 0) subInPlace(r, p);
    storeT(out, r, g);
#else
    uint carry = 0;
    UNROLL
    for (int i = T - 1; i >= 0; i--) {
        const size_t k = IDX(g, i);
        out[k] = addc(x[k], y[k], carry, &carry);
    }
    if (!carry) {
        int c = 0;
        for (int i = 0; i < T && c == 0; i++) {
            uint v = out[IDX(g, i)];
            if (v > p[i]) c =  1;
            if (v < p[i]) c = -1;
        }
        if (c < 0) return;
    }
    uint borrow = 0;
    for (int i = T - 1; i >= 0; i--) {
        const size_t k = IDX(g, i);
        out[k] = subb(out[k], p[i], borrow, &borrow);
    }
#endif
}

inline void op_submod(__global const uint *x, __global const uint *y,
                      __global uint *out, size_t g, const uint *p)
{
#if MPA_REGACC
    uint r[T];
    uint borrow = 0;
    UNROLL
    for (int i = T - 1; i >= 0; i--) {
        const size_t k = IDX(g, i);
        r[i] = subb(x[k], y[k], borrow, &borrow);
    }
    if (borrow) addInPlace(r, p);
    storeT(out, r, g);
#else
    uint borrow = 0;
    UNROLL
    for (int i = T - 1; i >= 0; i--) {
        const size_t k = IDX(g, i);
        out[k] = subb(x[k], y[k], borrow, &borrow);
    }
    if (borrow) {
        uint carry = 0;
        for (int i = T - 1; i >= 0; i--) {
            const size_t k = IDX(g, i);
            out[k] = addc(out[k], p[i], carry, &carry);
        }
    }
#endif
}

inline void op_mul_operand(__global const uint *x, __global const uint *y,
                           __global uint *out, size_t g)
{
#if MPA_REGACC
    uint acc[2 * T];
    UNROLL
    for (int i = 0; i < 2 * T; i++) acc[i] = 0;

    for (int i = T - 1; i >= 0; i--) {
        const uint xi = x[IDX(g, i)];
        uint U = 0;
        UNROLL
        for (int j = T - 1; j >= 0; j--) {
            uint hi, lo;
            mac(xi, y[IDX(g, j)], acc[i + j + 1], U, &hi, &lo);
            acc[i + j + 1] = lo;
            U = hi;
        }
        acc[i] = U;
    }
    UNROLL
    for (int i = 0; i < 2 * T; i++) out[IDX2(g, i)] = acc[i];
#else
    for (int i = 0; i < 2 * T; i++) out[IDX2(g, i)] = 0;
    for (int i = T - 1; i >= 0; i--) {
        const uint xi = x[IDX(g, i)];
        uint U = 0;
        for (int j = T - 1; j >= 0; j--) {
            uint hi, lo;
            mac(xi, y[IDX(g, j)], out[IDX2(g, i + j + 1)], U, &hi, &lo);
            out[IDX2(g, i + j + 1)] = lo;
            U = hi;
        }
        out[IDX2(g, i)] = U;
    }
#endif
}

inline void op_mul_product(__global const uint *x, __global const uint *y,
                           __global uint *out, size_t g)
{

#if MPA_MULHI
    uint c0 = 0, c1 = 0, c2 = 0;

    for (int k = 2 * T - 2; k >= 0; k--) {
        const int lo_i = (k - T + 1) > 0 ? (k - T + 1) : 0;
        const int hi_i = k < (T - 1) ? k : (T - 1);
        UNROLL
        for (int i = lo_i; i <= hi_i; i++) {
            uint hi, lo, carry;
            mac(x[IDX(g, i)], y[IDX(g, k - i)], 0, 0, &hi, &lo);
            c0 = addc(c0, lo, 0, &carry);
            c1 = addc(c1, hi, carry, &carry);
            c2 += carry;
        }
        out[IDX2(g, k + 1)] = c0;
        c0 = c1; c1 = c2; c2 = 0;
    }
    out[IDX2(g, 0)] = c0;
#else
    ulong accLo = 0;
    uint  accHi = 0;

    for (int k = 2 * T - 2; k >= 0; k--) {
        const int lo_i = (k - T + 1) > 0 ? (k - T + 1) : 0;
        const int hi_i = k < (T - 1) ? k : (T - 1);
        UNROLL
        for (int i = lo_i; i <= hi_i; i++) {
            const ulong prod = (ulong)x[IDX(g, i)] * (ulong)y[IDX(g, k - i)];
            accLo += prod;
            if (accLo < prod) accHi++;
        }
        out[IDX2(g, k + 1)] = (uint)accLo;
        accLo = (accLo >> 32) | (((ulong)accHi) << 32);
        accHi = 0;
    }
    out[IDX2(g, 0)] = (uint)accLo;
#endif
}


inline void copyN(uint *d, const uint *s)
{
    UNROLL
    for (int i = 0; i < T; i++) d[i] = s[i];
}

inline void zeroN(uint *a)
{
    UNROLL
    for (int i = 0; i < T; i++) a[i] = 0;
}

inline void oneN(uint *a)
{
    zeroN(a);
    a[T - 1] = 1;
}

inline int isZeroN(const uint *a)
{
    uint acc = 0;
    UNROLL
    for (int i = 0; i < T; i++) acc |= a[i];
    return acc == 0;
}

inline int cmpN(const uint *a, const uint *b)
{
    for (int i = 0; i < T; i++) {
        if (a[i] > b[i]) return  1;
        if (a[i] < b[i]) return -1;
    }
    return 0;
}

inline uint addN(uint *r, const uint *a, const uint *b)
{
    uint carry = 0;
    UNROLL
    for (int i = T - 1; i >= 0; i--) r[i] = addc(a[i], b[i], carry, &carry);
    return carry;
}

inline uint subN(uint *r, const uint *a, const uint *b)
{
    uint borrow = 0;
    UNROLL
    for (int i = T - 1; i >= 0; i--) r[i] = subb(a[i], b[i], borrow, &borrow);
    return borrow;
}

inline uint shl1N(uint *a)
{
    uint carry = 0;
    UNROLL
    for (int i = T - 1; i >= 0; i--) {
        const uint nx = a[i] >> 31;
        a[i] = (a[i] << 1) | carry;
        carry = nx;
    }
    return carry;
}

inline void shr1N(uint *a)
{
    uint carry = 0;
    UNROLL
    for (int i = 0; i < T; i++) {
        const uint nx = a[i] & 1u;
        a[i] = (a[i] >> 1) | (carry << 31);
        carry = nx;
    }
}

inline uint bitAtN(const uint *a, int bit)
{
    const int w = T - 1 - (bit >> 5);
    return (a[w] >> (bit & 31)) & 1u;
}

inline int bitLenN(const uint *a)
{
    for (int i = 0; i < T; i++) {
        if (a[i]) {
            int b = 0;
            uint v = a[i];
            while (v) { v >>= 1; b++; }
            return (T - 1 - i) * 32 + b;
        }
    }
    return 0;
}

inline void mulLowN(uint *r, const uint *a, const uint *b)
{
    uint acc[T];
    zeroN(acc);
    for (int i = T - 1; i >= 0; i--) {
        uint U = 0, hi, lo;
        for (int j = T - 1; j >= 0; j--) {
            const int k = i + j + 1 - T;
            if (k < 0) continue;
            mac(a[i], b[j], acc[k], U, &hi, &lo);
            acc[k] = lo;
            U = hi;
        }
    }
    copyN(r, acc);
}

inline void montMulPriv(uint *out, const uint *xbe, const uint *ybe,
                        const uint *nbe, uint m_prime)
{
    uint a[T], b[T], n[T], t[T + 2];

    UNROLL
    for (int i = 0; i < T; i++) {
        a[i] = xbe[T - 1 - i];
        b[i] = ybe[T - 1 - i];
        n[i] = nbe[T - 1 - i];
    }
    UNROLL
    for (int i = 0; i < T + 2; i++) t[i] = 0;

    for (int i = 0; i < T; i++) {
        uint C = 0, hi, lo, cc;
        const uint bi = b[i];
        UNROLL
        for (int j = 0; j < T; j++) {
            mac(a[j], bi, t[j], C, &hi, &lo);
            t[j] = lo; C = hi;
        }
        t[T]     = addc(t[T], C, 0, &cc);
        t[T + 1] = cc;

        const uint m = t[0] * m_prime;
        mac(m, n[0], t[0], 0, &hi, &lo);
        C = hi;
        UNROLL
        for (int j = 1; j < T; j++) {
            mac(m, n[j], t[j], C, &hi, &lo);
            t[j - 1] = lo; C = hi;
        }
        t[T - 1] = addc(t[T], C, 0, &cc);
        t[T]     = t[T + 1] + cc;
    }

    int ge = (t[T] != 0);
    if (!ge) {
        ge = 1;
        for (int j = T - 1; j >= 0; j--) {
            if (t[j] > n[j]) { ge = 1; break; }
            if (t[j] < n[j]) { ge = 0; break; }
        }
    }
    if (ge) {
        uint borrow = 0;
        UNROLL
        for (int j = 0; j < T; j++) t[j] = subb(t[j], n[j], borrow, &borrow);
    }

    UNROLL
    for (int i = 0; i < T; i++) out[i] = t[T - 1 - i];
}

inline void modDoubleN(uint *r, const uint *p)
{
    const uint c = shl1N(r);
    if (c || cmpN(r, p) >= 0) { uint tmp[T]; subN(tmp, r, p); copyN(r, tmp); }
}

inline void modAddN(uint *r, const uint *a, const uint *p)
{
    uint tmp[T];
    const uint c = addN(tmp, r, a);
    if (c || cmpN(tmp, p) >= 0) { uint t2[T]; subN(t2, tmp, p); copyN(r, t2); }
    else copyN(r, tmp);
}

inline void reduceN(uint *r, const uint *a, const uint *p)
{
    zeroN(r);
    for (int bit = 32 * T - 1; bit >= 0; bit--) {
        const uint c = shl1N(r);
        r[T - 1] |= bitAtN(a, bit);
        if (c || cmpN(r, p) >= 0) { uint tmp[T]; subN(tmp, r, p); copyN(r, tmp); }
    }
}

inline void divModN(uint *q, uint *r, const uint *a, const uint *b)
{
    zeroN(q);
    zeroN(r);
    if (isZeroN(b)) return;
    for (int bit = 32 * T - 1; bit >= 0; bit--) {
        const uint c = shl1N(r);
        r[T - 1] |= bitAtN(a, bit);
        shl1N(q);
        if (c || cmpN(r, b) >= 0) {
            uint tmp[T];
            subN(tmp, r, b);
            copyN(r, tmp);
            q[T - 1] |= 1u;
        }
    }
}

inline void toMontN(uint *out, const uint *a, const uint *p)
{
    uint acc[T];
    reduceN(acc, a, p);
    for (int i = 0; i < 32 * T; i++) modDoubleN(acc, p);
    copyN(out, acc);
}

inline void op_compare(__global const uint *x, __global const uint *y,
                       __global uint *out, size_t g)
{
    int c = 0;
    for (int i = 0; i < T; i++) {
        const uint xa = x[IDX(g, i)], yb = y[IDX(g, i)];
        if (xa > yb) { c =  1; break; }
        if (xa < yb) { c = -1; break; }
    }
    const uint fill = (c < 0) ? 0xFFFFFFFFu : 0u;
    UNROLL
    for (int i = 0; i < T; i++) out[IDX(g, i)] = fill;
    out[IDX(g, T - 1)] = (c < 0) ? 0xFFFFFFFFu : (uint)c;
}

inline void op_reduce(__global const uint *x, __global uint *out,
                      size_t g, const uint *p)
{
    uint a[T], r[T];
    UNROLL
    for (int i = 0; i < T; i++) a[i] = x[IDX(g, i)];
    reduceN(r, a, p);
    UNROLL
    for (int i = 0; i < T; i++) out[IDX(g, i)] = r[i];
}

inline void op_modmul(__global const uint *x, __global const uint *y,
                      __global uint *out, size_t g, const uint *p, uint m_prime)
{
    uint a[T], b[T], am[T], rm[T], r[T];
    UNROLL
    for (int i = 0; i < T; i++) { a[i] = x[IDX(g, i)]; b[i] = y[IDX(g, i)]; }
    toMontN(am, a, p);
    montMulPriv(rm, am, b, p, m_prime);
    copyN(r, rm);
    UNROLL
    for (int i = 0; i < T; i++) out[IDX(g, i)] = r[i];
}

inline void op_modmul_r2(__global const uint *x, __global const uint *y,
                         __global uint *out, size_t g, const uint *p,
                         __constant const uint *r2g, uint m_prime)
{
    uint a[T], b[T], r2[T], t1[T], r[T];
    UNROLL
    for (int i = 0; i < T; i++) {
        a[i]  = x[IDX(g, i)];
        b[i]  = y[IDX(g, i)];
        r2[i] = r2g[i];
    }
    montMulPriv(t1, a, b, p, m_prime);
    montMulPriv(r, t1, r2, p, m_prime);
    UNROLL
    for (int i = 0; i < T; i++) out[IDX(g, i)] = r[i];
}

inline void op_modexp(__global const uint *x, __global const uint *y,
                      __global uint *out, size_t g, const uint *p, uint m_prime)
{
    uint a[T], e[T], am[T], rm[T], one[T], r[T];
    UNROLL
    for (int i = 0; i < T; i++) { a[i] = x[IDX(g, i)]; e[i] = y[IDX(g, i)]; }

    oneN(one);
    toMontN(am, a, p);
    toMontN(rm, one, p);

    const int nb = bitLenN(e);
    for (int bit = nb - 1; bit >= 0; bit--) {
        uint t1[T];
        montMulPriv(t1, rm, rm, p, m_prime);
        copyN(rm, t1);
        if (bitAtN(e, bit)) {
            montMulPriv(t1, rm, am, p, m_prime);
            copyN(rm, t1);
        }
    }
    montMulPriv(r, rm, one, p, m_prime);
    UNROLL
    for (int i = 0; i < T; i++) out[IDX(g, i)] = r[i];
}

inline void op_exp(__global const uint *x, __global const uint *y,
                   __global uint *out, size_t g)
{
    uint a[T], e[T], r[T];
    UNROLL
    for (int i = 0; i < T; i++) { a[i] = x[IDX(g, i)]; e[i] = y[IDX(g, i)]; }
    oneN(r);
    const int nb = bitLenN(e);
    for (int bit = nb - 1; bit >= 0; bit--) {
        uint t1[T];
        mulLowN(t1, r, r);
        copyN(r, t1);
        if (bitAtN(e, bit)) {
            mulLowN(t1, r, a);
            copyN(r, t1);
        }
    }
    UNROLL
    for (int i = 0; i < T; i++) out[IDX(g, i)] = r[i];
}

inline void op_divide(__global const uint *x, __global const uint *y,
                      __global uint *out, size_t g)
{
    uint a[T], b[T], q[T], r[T];
    UNROLL
    for (int i = 0; i < T; i++) { a[i] = x[IDX(g, i)]; b[i] = y[IDX(g, i)]; }
    divModN(q, r, a, b);
    UNROLL
    for (int i = 0; i < T; i++) {
        out[IDX2(g, i)]     = q[i];
        out[IDX2(g, T + i)] = r[i];
    }
}

inline void op_isqrt(__global const uint *x, __global uint *out, size_t g)
{
    uint a[T], xc[T], yc[T], q[T], rr[T];
    UNROLL
    for (int i = 0; i < T; i++) a[i] = x[IDX(g, i)];

    if (isZeroN(a)) {
        UNROLL
        for (int i = 0; i < T; i++) out[IDX(g, i)] = 0;
        return;
    }

    const int nb = bitLenN(a);
    const int h = (nb + 1) >> 1;
    zeroN(xc);
    if (h >= 32 * T) copyN(xc, a);
    else xc[T - 1 - (h >> 5)] = 1u << (h & 31);

    for (int it = 0; it < 64 * T; it++) {
        divModN(q, rr, a, xc);
        addN(yc, xc, q);
        shr1N(yc);
        if (cmpN(yc, xc) >= 0) break;
        copyN(xc, yc);
    }
    UNROLL
    for (int i = 0; i < T; i++) out[IDX(g, i)] = xc[i];
}

inline void op_montgomery(__global const uint *x, __global const uint *y,
                          __global uint *out, size_t g,
                          const uint *pbe, uint m_prime)
{
#if MPA_FUSED_CIOS
    uint a[T], b[T], n[T], t[T + 2];

    UNROLL
    for (int i = 0; i < T; i++) {
        a[i] = x[IDX(g, T - 1 - i)];
        b[i] = y[IDX(g, T - 1 - i)];
        n[i] = pbe[T - 1 - i];
    }
    UNROLL
    for (int i = 0; i < T + 2; i++) t[i] = 0;

    for (int i = 0; i < T; i++) {
        uint C = 0, hi, lo, cc;
        const uint bi = b[i];

        UNROLL
        for (int j = 0; j < T; j++) {
            mac(a[j], bi, t[j], C, &hi, &lo);
            t[j] = lo; C = hi;
        }
        t[T]     = addc(t[T], C, 0, &cc);
        t[T + 1] = cc;

        const uint m = t[0] * m_prime;

        mac(m, n[0], t[0], 0, &hi, &lo);
        C = hi;
        UNROLL
        for (int j = 1; j < T; j++) {
            mac(m, n[j], t[j], C, &hi, &lo);
            t[j - 1] = lo; C = hi;
        }
        t[T - 1] = addc(t[T], C, 0, &cc);
        t[T]     = t[T + 1] + cc;
    }

    int ge = (t[T] != 0);
    if (!ge) {
        ge = 1;
        for (int j = T - 1; j >= 0; j--) {
            if (t[j] > n[j]) { ge = 1; break; }
            if (t[j] < n[j]) { ge = 0; break; }
        }
    }
    if (ge) {
        uint borrow = 0;
        UNROLL
        for (int j = 0; j < T; j++) t[j] = subb(t[j], n[j], borrow, &borrow);
    }

    UNROLL
    for (int i = 0; i < T; i++) out[IDX(g, i)] = t[T - 1 - i];

#else
    uint A[T + 1], xiy[T + 1], s[T + 2], mp[T + 1];

    UNROLL
    for (int i = 0; i <= T; i++) { A[i] = 0; xiy[i] = 0; s[i] = 0; }
    s[T + 1] = 0;

    const uint Yend = y[IDX(g, T - 1)];

    for (int i = T - 1; i >= 0; i--) {
        const uint xi = x[IDX(g, i)];
        const uint ui = (A[T] + xi * Yend) * m_prime;

        UNROLL
        for (int k = 0; k <= T; k++) xiy[k] = 0;
        {
            uint U = 0, hi, lo;
            UNROLL
            for (int k = T - 1; k >= 0; k--) {
                mac(y[IDX(g, k)], xi, xiy[k + 1], U, &hi, &lo);
                xiy[k + 1] = lo; U = hi;
            }
            xiy[0] = U;
        }
        {
            uint carry = 0;
            UNROLL
            for (int k = T; k >= 0; k--) s[k + 1] = addc(A[k], xiy[k], carry, &carry);
            s[0] = carry;
        }
        {
            uint U = 0, hi, lo;
            UNROLL
            for (int k = 0; k <= T; k++) mp[k] = 0;
            UNROLL
            for (int k = T - 1; k >= 0; k--) {
                mac(pbe[k], ui, mp[k + 1], U, &hi, &lo);
                mp[k + 1] = lo; U = hi;
            }
            mp[0] = U;
            uint carry = 0;
            UNROLL
            for (int k = T + 1; k >= 1; k--) s[k] = addc(s[k], mp[k - 1], carry, &carry);
            s[0] = s[0] + carry;
        }
        UNROLL
        for (int k = 0; k <= T; k++) A[k] = s[k];
    }

    int ge = (A[0] != 0);
    if (!ge) {
        ge = 1;
        for (int j = 0; j < T; j++) {
            if (A[j + 1] > pbe[j]) { ge = 1; break; }
            if (A[j + 1] < pbe[j]) { ge = 0; break; }
        }
    }
    if (ge) {
        uint borrow = 0;
        UNROLL
        for (int j = T - 1; j >= 0; j--) A[j + 1] = subb(A[j + 1], pbe[j], borrow, &borrow);
    }
    UNROLL
    for (int i = 0; i < T; i++) out[IDX(g, i)] = A[i + 1];
#endif
}

__kernel void mpaKernel(__global uint *input1, __global uint *input2,
                        __global uint *outputBytes,
                        __constant int *OPERATOR_WORDSIZE_BITSLENGHT_MPRIME,
                        __constant uint *globalPRIME)
{
    const size_t g = get_global_id(0);
    const uint m_prime = (uint)OPERATOR_WORDSIZE_BITSLENGHT_MPRIME[3];

    uint PRIME[T];
    UNROLL
    for (int i = 0; i < T; i++) PRIME[i] = globalPRIME[i];

    switch (OPERATOR_WORDSIZE_BITSLENGHT_MPRIME[0]) {
    case ADD:
        op_add(input1, input2, outputBytes, g); break;
    case SUBTRACT:
        op_sub(input1, input2, outputBytes, g); break;
    case ADDMOD:
        op_addmod(input1, input2, outputBytes, g, PRIME); break;
    case SUBTRACTMOD:
        op_submod(input1, input2, outputBytes, g, PRIME); break;
    case MULTIPLYOPERANDSCANNING:
        op_mul_operand(input1, input2, outputBytes, g); break;
    case MULTIPLYPRODUCTSCANNING:
        op_mul_product(input1, input2, outputBytes, g); break;
    case MONTGOMERYMULTIPLICATION:
        op_montgomery(input1, input2, outputBytes, g, PRIME, m_prime); break;
    case COMPARE:
        op_compare(input1, input2, outputBytes, g); break;
    case REDUCE:
        op_reduce(input1, outputBytes, g, PRIME); break;
    case MODMUL:
        op_modmul(input1, input2, outputBytes, g, PRIME, m_prime); break;
    case MODMUL_R2:
        op_modmul_r2(input1, input2, outputBytes, g, PRIME,
                     globalPRIME + WORDLENGTH_T, m_prime); break;
    case MODEXP:
        op_modexp(input1, input2, outputBytes, g, PRIME, m_prime); break;
    case EXPONENTIATION:
        op_exp(input1, input2, outputBytes, g); break;
    case DIVIDE:
        op_divide(input1, input2, outputBytes, g); break;
    case ISQRT:
        op_isqrt(input1, outputBytes, g); break;
    default:
        for (int i = 0; i < T; i++) outputBytes[IDX(g, i)] = 0xFFFFFFFFu;
        break;
    }
}
