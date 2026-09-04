#define T 8

#if MPA_UNROLL
#define UNROLL __attribute__((opencl_unroll_hint))
#else
#define UNROLL
#endif

inline uint addc(uint a, uint b, uint cin, uint *cout)
{
    ulong s = (ulong)a + (ulong)b + (ulong)cin;
    *cout = (uint)(s >> 32);
    return (uint)s;
}

inline uint subb(uint a, uint b, uint bin, uint *bout)
{
    ulong d = (ulong)a - (ulong)b - (ulong)bin;
    *bout = (uint)((d >> 32) & 1u);
    return (uint)d;
}

inline void mac(uint a, uint b, uint c, uint d, uint *hi, uint *lo)
{
    ulong p = (ulong)a * (ulong)b + (ulong)c + (ulong)d;
    *hi = (uint)(p >> 32);
    *lo = (uint)p;
}

inline void cpy8(uint *d, const uint *s)
{
    UNROLL
    for (int i = 0; i < T; i++) d[i] = s[i];
}

inline void zero8(uint *d)
{
    UNROLL
    for (int i = 0; i < T; i++) d[i] = 0;
}

inline int isz8(const uint *a)
{
    uint acc = 0;
    UNROLL
    for (int i = 0; i < T; i++) acc |= a[i];
    return acc == 0;
}

inline int cmp8c(const uint *a, __constant const uint *b)
{
    for (int i = 0; i < T; i++) {
        if (a[i] > b[i]) return  1;
        if (a[i] < b[i]) return -1;
    }
    return 0;
}

inline int eq8(const uint *a, const uint *b)
{
    uint d = 0;
    UNROLL
    for (int i = 0; i < T; i++) d |= a[i] ^ b[i];
    return d == 0;
}

inline void modaddC(uint *r, const uint *a, const uint *b, __constant const uint *m)
{
    uint c = 0, t[T];
    UNROLL
    for (int i = T - 1; i >= 0; i--) t[i] = addc(a[i], b[i], c, &c);
    if (c || cmp8c(t, m) >= 0) {
        uint br = 0;
        UNROLL
        for (int i = T - 1; i >= 0; i--) t[i] = subb(t[i], m[i], br, &br);
    }
    cpy8(r, t);
}

inline void modsubC(uint *r, const uint *a, const uint *b, __constant const uint *m)
{
    uint br = 0, t[T];
    UNROLL
    for (int i = T - 1; i >= 0; i--) t[i] = subb(a[i], b[i], br, &br);
    if (br) {
        uint c = 0;
        UNROLL
        for (int i = T - 1; i >= 0; i--) t[i] = addc(t[i], m[i], c, &c);
    }
    cpy8(r, t);
}

inline void montmulC(uint *out, const uint *xbe, const uint *ybe,
                     __constant const uint *nbe, uint mprime)
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

        const uint mm = t[0] * mprime;
        mac(mm, n[0], t[0], 0, &hi, &lo);
        C = hi;
        UNROLL
        for (int j = 1; j < T; j++) {
            mac(mm, n[j], t[j], C, &hi, &lo);
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
        uint br = 0;
        UNROLL
        for (int j = 0; j < T; j++) t[j] = subb(t[j], n[j], br, &br);
    }
    UNROLL
    for (int i = 0; i < T; i++) out[i] = t[T - 1 - i];
}

inline void montmulCC(uint *out, const uint *xbe, __constant const uint *ybe,
                      __constant const uint *nbe, uint mprime)
{
    uint yy[T];
    UNROLL
    for (int i = 0; i < T; i++) yy[i] = ybe[i];
    montmulC(out, xbe, yy, nbe, mprime);
}

inline uint bit8(const uint *a, int b)
{
    return (a[T - 1 - (b >> 5)] >> (b & 31)) & 1u;
}

typedef struct { uint X[T], Y[T], Z[T]; } Jac;

inline int jIsInf(const Jac *P) { return isz8(P->Z); }

inline void jDouble(Jac *R, const Jac *P, __constant const uint *p, uint mp)
{
    if (jIsInf(P) || isz8(P->Y)) { zero8(R->Z); return; }

    uint delta[T], gamma[T], beta[T], alpha[T], t1[T], t2[T], t3[T], X3[T];

    montmulC(delta, P->Z, P->Z, p, mp);
    montmulC(gamma, P->Y, P->Y, p, mp);
    montmulC(beta,  P->X, gamma, p, mp);

    modsubC(t1, P->X, delta, p);
    modaddC(t2, P->X, delta, p);
    montmulC(t3, t1, t2, p, mp);
    modaddC(alpha, t3, t3, p);
    modaddC(alpha, alpha, t3, p);

    montmulC(t1, alpha, alpha, p, mp);
    modaddC(t2, beta, beta, p);
    modaddC(t2, t2, t2, p);
    modaddC(t3, t2, t2, p);
    modsubC(X3, t1, t3, p);

    modaddC(t3, P->Y, P->Z, p);
    montmulC(t3, t3, t3, p, mp);
    modsubC(t3, t3, gamma, p);
    modsubC(t3, t3, delta, p);

    modsubC(t1, t2, X3, p);
    montmulC(t1, alpha, t1, p, mp);
    montmulC(t2, gamma, gamma, p, mp);
    modaddC(t2, t2, t2, p);
    modaddC(t2, t2, t2, p);
    modaddC(t2, t2, t2, p);
    modsubC(t1, t1, t2, p);

    cpy8(R->X, X3);
    cpy8(R->Y, t1);
    cpy8(R->Z, t3);
}

inline void jAddAffine(Jac *R, const Jac *P, const uint *Qx, const uint *Qy,
                       __constant const uint *p, uint mp, __constant const uint *r2p)
{
    if (jIsInf(P)) {
        cpy8(R->X, Qx);
        cpy8(R->Y, Qy);
        zero8(R->Z); R->Z[T - 1] = 1;
        montmulCC(R->Z, R->Z, r2p, p, mp);
        return;
    }

    uint Z1Z1[T], U2[T], S2[T], H[T], HH[T], I[T], J[T], rr[T], V[T], t1[T], X3[T], Y3[T];

    montmulC(Z1Z1, P->Z, P->Z, p, mp);
    montmulC(U2, Qx, Z1Z1, p, mp);
    montmulC(t1, P->Z, Z1Z1, p, mp);
    montmulC(S2, Qy, t1, p, mp);

    modsubC(H, U2, P->X, p);
    modsubC(rr, S2, P->Y, p);

    if (isz8(H)) {
        if (isz8(rr)) { jDouble(R, P, p, mp); return; }
        zero8(R->Z);
        return;
    }

    montmulC(HH, H, H, p, mp);
    modaddC(I, HH, HH, p);
    modaddC(I, I, I, p);
    montmulC(J, H, I, p, mp);
    modaddC(rr, rr, rr, p);
    montmulC(V, P->X, I, p, mp);

    montmulC(X3, rr, rr, p, mp);
    modsubC(X3, X3, J, p);
    modsubC(X3, X3, V, p);
    modsubC(X3, X3, V, p);

    modsubC(Y3, V, X3, p);
    montmulC(Y3, rr, Y3, p, mp);
    montmulC(t1, P->Y, J, p, mp);
    modaddC(t1, t1, t1, p);
    modsubC(Y3, Y3, t1, p);

    modaddC(t1, P->Z, H, p);
    montmulC(t1, t1, t1, p, mp);
    modsubC(t1, t1, Z1Z1, p);
    modsubC(t1, t1, HH, p);

    cpy8(R->X, X3);
    cpy8(R->Y, Y3);
    cpy8(R->Z, t1);
}

inline void jAddAffineC(Jac *R, const Jac *P, __constant const uint *Qx,
                        __constant const uint *Qy, __constant const uint *p,
                        uint mp, __constant const uint *r2p)
{
    uint qx[T], qy[T];
    UNROLL
    for (int i = 0; i < T; i++) { qx[i] = Qx[i]; qy[i] = Qy[i]; }
    jAddAffine(R, P, qx, qy, p, mp, r2p);
}

inline void modexpMontC(uint *out, const uint *baseM, const uint *e,
                        __constant const uint *m, __constant const uint *r2, uint mp)
{
    uint acc[T], t[T], one[T];
    zero8(one); one[T - 1] = 1;
    montmulCC(acc, one, r2, m, mp);

    int top = -1;
    for (int i = 32 * T - 1; i >= 0; i--) if (bit8(e, i)) { top = i; break; }

    for (int i = top; i >= 0; i--) {
        montmulC(t, acc, acc, m, mp);
        cpy8(acc, t);
        if (bit8(e, i)) {
            montmulC(t, acc, baseM, m, mp);
            cpy8(acc, t);
        }
    }
    cpy8(out, acc);
}

__kernel void ecdsaVerify(__global const uint *sig,
                          __global const uint *msg,
                          __global const uint *pub,
                          __global uchar *valid,
                          __constant uint *par)
{
    const size_t g = get_global_id(0);

    __constant const uint *p    = par;
    __constant const uint *r2p  = par + T;
    __constant const uint *n    = par + 2 * T;
    __constant const uint *r2n  = par + 3 * T;
    __constant const uint *Gx   = par + 4 * T;
    __constant const uint *Gy   = par + 5 * T;
    __constant const uint *nmod = par + 8 * T;
    const uint mp = par[6 * T];
    const uint mn = par[6 * T + 1];
    const uint rPlusNFits = par[6 * T + 2];

    uint r[T], s[T], e[T], Qx[T], Qy[T];
    UNROLL
    for (int i = 0; i < T; i++) {
        r[i]  = sig[g * 2 * T + i];
        s[i]  = sig[g * 2 * T + T + i];
        e[i]  = msg[g * T + i];
        Qx[i] = pub[g * 2 * T + i];
        Qy[i] = pub[g * 2 * T + T + i];
    }

    valid[g] = 0;

    if (isz8(r) || isz8(s)) return;
    if (cmp8c(r, n) >= 0 || cmp8c(s, n) >= 0) return;

    if (cmp8c(e, n) >= 0) {
        uint br = 0;
        UNROLL
        for (int i = T - 1; i >= 0; i--) e[i] = subb(e[i], n[i], br, &br);
    }

    uint u1[T], u2[T];
    {
        uint sM[T], w[T], nm2[T], one[T];
        zero8(one); one[T - 1] = 1;
        montmulCC(sM, s, r2n, n, mn);
        {
            uint two[T], br = 0;
            zero8(two); two[T - 1] = 2;
            UNROLL
            for (int i = 0; i < T; i++) nm2[i] = n[i];
            UNROLL
            for (int i = T - 1; i >= 0; i--) nm2[i] = subb(nm2[i], two[i], br, &br);
        }
        modexpMontC(w, sM, nm2, n, r2n, mn);

        montmulCC(u1, e, r2n, n, mn);
        montmulC(u1, u1, w, n, mn);
        montmulC(u1, u1, one, n, mn);

        montmulCC(u2, r, r2n, n, mn);
        montmulC(u2, u2, w, n, mn);
        montmulC(u2, u2, one, n, mn);
    }

    montmulCC(Qx, Qx, r2p, p, mp);
    montmulCC(Qy, Qy, r2p, p, mp);

    Jac acc;
    zero8(acc.X); zero8(acc.Y); zero8(acc.Z);

    for (int i = 32 * T - 1; i >= 0; i--) {
        jDouble(&acc, &acc, p, mp);
        if (bit8(u1, i)) jAddAffineC(&acc, &acc, Gx, Gy, p, mp, r2p);
        if (bit8(u2, i)) jAddAffine(&acc, &acc, Qx, Qy, p, mp, r2p);
    }

    if (jIsInf(&acc)) return;

    {
        uint z2[T], t[T], rM[T];
        montmulC(z2, acc.Z, acc.Z, p, mp);

        montmulCC(rM, r, r2p, p, mp);
        montmulC(t, rM, z2, p, mp);
        if (eq8(t, acc.X)) { valid[g] = 1; return; }

        if (rPlusNFits) {
            uint rn[T];
            uint c = 0;
            UNROLL
            for (int i = T - 1; i >= 0; i--) rn[i] = addc(r[i], nmod[i], c, &c);
            if (!c && cmp8c(rn, p) < 0) {
                montmulCC(rM, rn, r2p, p, mp);
                montmulC(t, rM, z2, p, mp);
                if (eq8(t, acc.X)) valid[g] = 1;
            }
        }
    }
}
