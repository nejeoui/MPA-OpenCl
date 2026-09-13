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

inline int cmp8(const uint *a, const uint *b)
{
    for (int i = 0; i < T; i++) {
        if (a[i] > b[i]) return  1;
        if (a[i] < b[i]) return -1;
    }
    return 0;
}

inline void modadd(uint *r, const uint *a, const uint *b, const uint *m)
{
    uint c = 0, t[T];
    UNROLL
    for (int i = T - 1; i >= 0; i--) t[i] = addc(a[i], b[i], c, &c);
    if (c || cmp8(t, m) >= 0) {
        uint br = 0;
        UNROLL
        for (int i = T - 1; i >= 0; i--) t[i] = subb(t[i], m[i], br, &br);
    }
    cpy8(r, t);
}

inline void modsub(uint *r, const uint *a, const uint *b, const uint *m)
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

inline void montmul(uint *out, const uint *xbe, const uint *ybe,
                    const uint *nbe, uint mprime)
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

inline uint bit8(const uint *a, int b)
{
    return (a[T - 1 - (b >> 5)] >> (b & 31)) & 1u;
}

typedef struct { uint X[T], Y[T], Z[T]; } Jac;

inline int jIsInf(const Jac *P) { return isz8(P->Z); }

inline void jSetInf(Jac *P)
{
    zero8(P->X); zero8(P->Y); zero8(P->Z);
    P->X[T - 1] = 1;
    P->Y[T - 1] = 1;
}

inline void jDouble(Jac *R, const Jac *P, const uint *p, uint mp)
{
    if (jIsInf(P) || isz8(P->Y)) { jSetInf(R); return; }

    uint delta[T], gamma[T], beta[T], alpha[T], t1[T], t2[T], t3[T];

    montmul(delta, P->Z, P->Z, p, mp);
    montmul(gamma, P->Y, P->Y, p, mp);
    montmul(beta,  P->X, gamma, p, mp);

    modsub(t1, P->X, delta, p);
    modadd(t2, P->X, delta, p);
    montmul(t3, t1, t2, p, mp);
    modadd(alpha, t3, t3, p);
    modadd(alpha, alpha, t3, p);

    montmul(t1, alpha, alpha, p, mp);
    modadd(t2, beta, beta, p);
    modadd(t2, t2, t2, p);
    modadd(t3, t2, t2, p);
    modsub(t1, t1, t3, p);

    modadd(t3, P->Y, P->Z, p);
    montmul(t3, t3, t3, p, mp);
    modsub(t3, t3, gamma, p);
    modsub(t3, t3, delta, p);

    uint X3[T];
    cpy8(X3, t1);

    modsub(t1, t2, X3, p);
    montmul(t1, alpha, t1, p, mp);
    montmul(t2, gamma, gamma, p, mp);
    modadd(t2, t2, t2, p);
    modadd(t2, t2, t2, p);
    modadd(t2, t2, t2, p);
    modsub(t1, t1, t2, p);

    cpy8(R->X, X3);
    cpy8(R->Y, t1);
    cpy8(R->Z, t3);
}

inline void jAdd(Jac *R, const Jac *P, const Jac *Q, const uint *p, uint mp)
{
    if (jIsInf(P)) { cpy8(R->X, Q->X); cpy8(R->Y, Q->Y); cpy8(R->Z, Q->Z); return; }
    if (jIsInf(Q)) { cpy8(R->X, P->X); cpy8(R->Y, P->Y); cpy8(R->Z, P->Z); return; }

    uint Z1Z1[T], Z2Z2[T], U1[T], U2[T], S1[T], S2[T], H[T], I[T], J[T], rr[T], V[T], t1[T];

    montmul(Z1Z1, P->Z, P->Z, p, mp);
    montmul(Z2Z2, Q->Z, Q->Z, p, mp);
    montmul(U1, P->X, Z2Z2, p, mp);
    montmul(U2, Q->X, Z1Z1, p, mp);
    montmul(t1, Q->Z, Z2Z2, p, mp);
    montmul(S1, P->Y, t1, p, mp);
    montmul(t1, P->Z, Z1Z1, p, mp);
    montmul(S2, Q->Y, t1, p, mp);

    modsub(H, U2, U1, p);
    modsub(rr, S2, S1, p);

    if (isz8(H)) {
        if (isz8(rr)) { jDouble(R, P, p, mp); return; }
        jSetInf(R);
        return;
    }

    modadd(I, H, H, p);
    montmul(I, I, I, p, mp);
    montmul(J, H, I, p, mp);
    modadd(rr, rr, rr, p);
    montmul(V, U1, I, p, mp);

    uint X3[T], Y3[T], Z3[T];
    montmul(X3, rr, rr, p, mp);
    modsub(X3, X3, J, p);
    modsub(X3, X3, V, p);
    modsub(X3, X3, V, p);

    modsub(Y3, V, X3, p);
    montmul(Y3, rr, Y3, p, mp);
    montmul(t1, S1, J, p, mp);
    modadd(t1, t1, t1, p);
    modsub(Y3, Y3, t1, p);

    modadd(Z3, P->Z, Q->Z, p);
    montmul(Z3, Z3, Z3, p, mp);
    modsub(Z3, Z3, Z1Z1, p);
    modsub(Z3, Z3, Z2Z2, p);
    montmul(Z3, Z3, H, p, mp);

    cpy8(R->X, X3);
    cpy8(R->Y, Y3);
    cpy8(R->Z, Z3);
}

inline void modexpMont(uint *out, const uint *baseM, const uint *e,
                       const uint *m, const uint *r2, uint mp)
{
    uint acc[T], t[T], one[T];
    zero8(one); one[T - 1] = 1;
    montmul(acc, one, r2, m, mp);

    int top = -1;
    for (int i = 32 * T - 1; i >= 0; i--) if (bit8(e, i)) { top = i; break; }

    for (int i = top; i >= 0; i--) {
        montmul(t, acc, acc, m, mp);
        cpy8(acc, t);
        if (bit8(e, i)) {
            montmul(t, acc, baseM, m, mp);
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

    uint p[T], r2p[T], n[T], r2n[T], Gx[T], Gy[T];
    UNROLL
    for (int i = 0; i < T; i++) {
        p[i]   = par[i];
        r2p[i] = par[T + i];
        n[i]   = par[2 * T + i];
        r2n[i] = par[3 * T + i];
        Gx[i]  = par[4 * T + i];
        Gy[i]  = par[5 * T + i];
    }
    const uint mp = par[6 * T];
    const uint mn = par[6 * T + 1];

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
    if (cmp8(r, n) >= 0 || cmp8(s, n) >= 0) return;

    if (cmp8(e, n) >= 0) {
        uint br = 0;
        UNROLL
        for (int i = T - 1; i >= 0; i--) e[i] = subb(e[i], n[i], br, &br);
    }

    uint sM[T], w[T], nm2[T], one[T];
    zero8(one); one[T - 1] = 1;

    montmul(sM, s, r2n, n, mn);

    {
        uint two[T], br = 0;
        zero8(two); two[T - 1] = 2;
        cpy8(nm2, n);
        UNROLL
        for (int i = T - 1; i >= 0; i--) nm2[i] = subb(nm2[i], two[i], br, &br);
    }

    modexpMont(w, sM, nm2, n, r2n, mn);

    uint eM[T], rM[T], u1[T], u2[T];
    montmul(eM, e, r2n, n, mn);
    montmul(rM, r, r2n, n, mn);
    montmul(u1, eM, w, n, mn);
    montmul(u2, rM, w, n, mn);
    montmul(u1, u1, one, n, mn);
    montmul(u2, u2, one, n, mn);

    uint QxM[T], QyM[T];
    montmul(QxM, Qx, r2p, p, mp);
    montmul(QyM, Qy, r2p, p, mp);

    Jac Pg, Pq, Pgq, acc, tmp;
    zero8(Pg.Z); Pg.Z[T - 1] = 1;
    montmul(Pg.Z, Pg.Z, r2p, p, mp);
    cpy8(Pg.X, Gx);
    cpy8(Pg.Y, Gy);

    cpy8(Pq.X, QxM);
    cpy8(Pq.Y, QyM);
    cpy8(Pq.Z, Pg.Z);

    jAdd(&Pgq, &Pg, &Pq, p, mp);

    zero8(acc.X); acc.X[T - 1] = 1;
    montmul(acc.X, acc.X, r2p, p, mp);
    cpy8(acc.Y, acc.X);
    zero8(acc.Z);

    for (int i = 32 * T - 1; i >= 0; i--) {
        jDouble(&tmp, &acc, p, mp);
        acc = tmp;
        const uint b1 = bit8(u1, i), b2 = bit8(u2, i);
        if (b1 && b2)       { jAdd(&tmp, &acc, &Pgq, p, mp); acc = tmp; }
        else if (b1)        { jAdd(&tmp, &acc, &Pg,  p, mp); acc = tmp; }
        else if (b2)        { jAdd(&tmp, &acc, &Pq,  p, mp); acc = tmp; }
    }

    if (jIsInf(&acc)) return;

    uint zinv[T], pm2[T], z2[T], xa[T];
    {
        uint two[T];
        zero8(two); two[T - 1] = 2;
        cpy8(pm2, p);
        uint b2 = 0;
        UNROLL
        for (int i = T - 1; i >= 0; i--) pm2[i] = subb(pm2[i], two[i], b2, &b2);
    }
    modexpMont(zinv, acc.Z, pm2, p, r2p, mp);
    montmul(z2, zinv, zinv, p, mp);
    montmul(xa, acc.X, z2, p, mp);
    montmul(xa, xa, one, p, mp);

    if (cmp8(xa, n) >= 0) {
        uint b3 = 0;
        UNROLL
        for (int i = T - 1; i >= 0; i--) xa[i] = subb(xa[i], n[i], b3, &b3);
    }

    valid[g] = (cmp8(xa, r) == 0) ? 1 : 0;
}
