#define TWOPOW_W 256
#define ADD 1
#define SUBTRACT 2
#define ADDMOD 3
#define SUBTRACTMOD 4
#define MULTIPLYOPERANDSCANNING 5
#define MULTIPLYPRODUCTSCANNING 6
#define MONTGOMERYMULTIPLICATION 7
#define COMPARE        8
#define REDUCE         9
#define MODMUL        10
#define MODEXP        11
#define EXPONENTIATION 12
#define DIVIDE        13
#define ISQRT         14
#define MODMUL_R2     15

typedef unsigned char xword;
#define XW      8
#define XWMASK  0xFFu
#define XT      WORDLENGTH_T
#define XIDX(g,i)   ((size_t)(g) * (size_t)XT + (size_t)(i))
#define XIDX2(g,i)  ((size_t)(g) * (size_t)(2 * XT) + (size_t)(i))

#include <mpaKernels_8bits.h>

void addPrime(__global unsigned char*  outputBytes, const size_t ID, __private unsigned char PRIME[]){
    uint carry = 0;
    int i=0;
    for ( i = WORDLENGTH_T-1; i >= 0; i--)
    {
        const size_t index=ID*WORDLENGTH_T+i;
        const uint somme = (uint)outputBytes[index] + (uint)PRIME[i] + (uint)carry;
        outputBytes[index] = (unsigned char)somme;
        carry = (uint)(somme >> 8);
    }
}

void subtractPrime(__global unsigned char*  outputBytes, const size_t ID,__private unsigned char PRIME[]){
    uint borrow = 0;
    int i=0;
    for ( i = WORDLENGTH_T-1; i >= 0; i--) {
        const size_t index=ID*WORDLENGTH_T+i;
        const uint diff = (uint)outputBytes[index] - (uint)PRIME[i] - (uint)borrow;
        outputBytes[index] = (unsigned char)diff;
        borrow = (uint)((diff >> 8) & 1u);
    }
}

char compareWithPrime(__global unsigned char*  outputBytes, const size_t ID, __private unsigned char PRIME[]){

int i=0;
for ( i = 0; i < WORDLENGTH_T; i++)
{
    const size_t index=ID*WORDLENGTH_T+i;
        if (outputBytes[index] > PRIME[i])
            return 1;
        if (outputBytes[index] < PRIME[i])
            return -1;
    }
    return 0;

}

void add(__global unsigned char* input1, __global unsigned char* input2, __global unsigned char* outputBytes, const size_t ID)
{
    uint carry = 0;
    int i=0;
    for (i=WORDLENGTH_T-1; i >= 0 ; i--)
    {
        const size_t index=ID*WORDLENGTH_T+i;
        const uint somme = (uint)carry + (uint)input1[index] + (uint)input2[index];
        outputBytes[index] = (unsigned char)somme;
        carry = (uint)(somme >> 8);
    }
}

void multiplyOperandScanning(__global unsigned char* input1, __global unsigned char* input2, __global unsigned char* outputBytes, const size_t ID)
{
    unsigned long int UV=0;
    unsigned long int U=0;
    unsigned long int V=0;
    int i;
    for( i=2*WORDLENGTH_T-1;i>=0;i--)
        outputBytes[ID*2*WORDLENGTH_T+i]=0;

    for( i=WORDLENGTH_T-1;i>=0;i--) {
        U=0;
        const size_t indexI=ID*WORDLENGTH_T+i;
        int j;
        for( j=WORDLENGTH_T-1;j>=0;j--) {
            const size_t indexJ=ID*WORDLENGTH_T+j;
            UV=(uint)outputBytes[ID*2*WORDLENGTH_T+i+j+1]+((uint)input1[indexI])*((uint)input2[indexJ])+U;
            U=(UV&0xFF00)>>8;
            V=UV&0xFF;

            outputBytes[ID*2*WORDLENGTH_T+ i+j+1]=(unsigned char)V;
        }
        outputBytes[ID*2*WORDLENGTH_T+ i]=(unsigned char)U;
    }
}
int MIN(int x,int y) {
    if(x<y) return x;
    else return y;
}
 int MAX(int x,int y) {
    if(x>y) return x;
    else return y;
}
void multiplyProductScanning(__global unsigned char* input1, __global unsigned char* input2, __global unsigned char* outputBytes, const size_t ID)
{
    ulong carry = 0;
    int k;
   for( k=2*WORDLENGTH_T-2;k>=0;k--) {
        ulong acc = carry;
        int i;
        for( i=MAX(0,k-WORDLENGTH_T+1);i<=MIN(k,WORDLENGTH_T-1);i++)
            acc += ((ulong)input1[ID*WORDLENGTH_T+i]) * ((ulong)input2[ID*WORDLENGTH_T+k-i]);
        outputBytes[ID*2*WORDLENGTH_T+ k+1]=(unsigned char)(acc & 0xFF);
        carry = acc >> 8;
    }
    outputBytes[ID*2*WORDLENGTH_T]=(unsigned char)carry;
}

void subtractPositive(__global unsigned char* input1, __global unsigned char* input2, __global unsigned char* outputBytes, const size_t ID)
{
    uint borrow = 0;
    int i=0;
    for (i = WORDLENGTH_T-1; i >= 0; i--) {
         const size_t index=ID*WORDLENGTH_T+i;
        const uint diff = (uint)input1[index] - (uint)input2[index] - (uint)borrow;
        outputBytes[index] = (unsigned char)diff;
        borrow = (uint)((diff >> 8) & 1u);
    }
}

 void addMod(__global unsigned char* input1, __global unsigned char* input2, __global unsigned char* outputBytes, const size_t ID, __private unsigned char PRIME[])
{
    uint carry = 0;
    int i=0;

    for ( i = WORDLENGTH_T-1; i >= 0; i--)
    {
       const size_t index=ID*WORDLENGTH_T+i;
        const uint somme = (uint)input1[index] + (uint)input2[index] + (uint)carry;
        outputBytes[index] = (unsigned char)somme;
        carry = (uint)(somme >> 8);
    }

    if (carry == 1) {
        subtractPrime(outputBytes,ID,PRIME);
    }
    else if(compareWithPrime( outputBytes,ID,PRIME)>=0) {
        subtractPrime(outputBytes,ID,PRIME);
    }
}

 void subtractMod(__global unsigned char* input1, __global unsigned char* input2, __global unsigned char* outputBytes, const size_t ID,  __private unsigned char PRIME[])
{
    uint borrow = 0;
    int i=0;
    for ( i = WORDLENGTH_T-1; i >= 0; i--)
    {
        const size_t  index=ID*WORDLENGTH_T+i;
        const uint diff = (uint)input1[index] - (uint)input2[index] - (uint)borrow;
        outputBytes[index] = (unsigned char)diff;
        borrow = (uint)((diff >> 8) & 1u);
    }
    if (borrow == 1) {
        addPrime(outputBytes,ID,PRIME);
    }
}

  void montgomeryMultiplication(__global unsigned char*  x,__global unsigned char* y,__global unsigned char* result,const size_t ID,__private unsigned char PRIME[],const unsigned int m_prime) {
    __private   unsigned char resultPrivate[WORDLENGTH_T+1];
    __private   unsigned char xiy[WORDLENGTH_T+1];
    __private   unsigned char Aplusxiy[WORDLENGTH_T+2];

    __private   unsigned char cteUI[WORDLENGTH_T+1];
    int i;
    const unsigned char Yend=y[ID*WORDLENGTH_T+WORDLENGTH_T-1];
     for( i=WORDLENGTH_T;i>=0;i--) {
resultPrivate[i]=0;
            xiy[i]=0;
       Aplusxiy[i]=0;

 }

 Aplusxiy[WORDLENGTH_T+1]=0;
    for( i=WORDLENGTH_T-1;i>=0;i--) {
        size_t xindex=i+WORDLENGTH_T*ID;
        unsigned int ui=(((uint)resultPrivate[WORDLENGTH_T]+((uint)x[xindex])*Yend)*m_prime)&0xFF;
        multiplyNoOverFlow1xWORDLENGTH(x[xindex],ID,y,xiy);
        addNoOverFlowPrivate_XIY(resultPrivate,xiy,Aplusxiy);
        addNoOverFlowPrivateAplusxiy(ui,Aplusxiy,cteUI,PRIME);
        rightShiftFormby1InResultPriv(Aplusxiy,resultPrivate);

    }
    if(compareResultPrivPrime(resultPrivate,PRIME)>=0) subtractPositiveResultPrivate(resultPrivate,PRIME);
    copyResultPrivTo(result,resultPrivate,ID);
}

 void rightShiftFormby1InResultPriv(__private  unsigned char Aplusxiy[],__private  unsigned char resultPrivate[]) {
    int i;
    for( i=0;i<WORDLENGTH_T+1;i++) resultPrivate[i]=Aplusxiy[i];
}
void subtractPositiveResultPrivate(__private unsigned char resultPrivate[],__private unsigned char PRIME[]){
    uint borrow = 0;
    int i;
    for ( i = WORDLENGTH_T-1; i >= 0; i--) {
        const uint diff = (uint)resultPrivate[i+1] - (uint)PRIME[i] - (uint)borrow;
        resultPrivate[i+1] = (unsigned char)diff;
        borrow = (uint)((diff >> 8) & 1u);
    }
    resultPrivate[0]=(unsigned char)(resultPrivate[0]-borrow);
}
int compareResultPrivPrime(__private unsigned char resultPrivate[],__private unsigned char PRIME[]){
    int i;
    if (resultPrivate[0]!=0) return 1;
    else
    for ( i = 0; i < WORDLENGTH_T; i++) {
        if (resultPrivate[i+1] >PRIME[i])
            return 1;
        if (resultPrivate[i+1] < PRIME[i])
            return -1;
    }
    return 0;
}

 void multiplyNoOverFlow1xWORDLENGTH(unsigned char n,const size_t ID,__global unsigned char* y,__private unsigned char xiy[]) {

        int alength=WORDLENGTH_T;
    unsigned long int UV=0;
    unsigned long int U=0;
    unsigned long int V=0;
    int i;
        for( i=alength;i>=0;i--)
            xiy[i]=0;

        for( i=alength-1;i>=0;i--)
        {
            U=0;
             {
                UV=(uint)xiy[i+1]+ ((uint)y[WORDLENGTH_T*ID+i])*n + U;
                U=(UV&0xFF00)>>8;
                V=UV&0xFF;

                xiy[i+1]=(unsigned char)V;
            }
            xiy[i]=(unsigned char)U;
        }

}

 void addNoOverFlowPrivate_XIY(__private unsigned char resultPrivate[],__private unsigned char xiy[],__private unsigned char Aplusxiy[]) {
   unsigned int carry = 0;
   unsigned  int somme = 0;

    int i;
    for ( i = WORDLENGTH_T; i >= 0; i--)
    {
        somme = (uint)resultPrivate[i] + (uint)xiy[i]  + carry;
        carry = 0;
        if (somme >= TWOPOW_W) {
            somme -= TWOPOW_W;
            carry = 1;
        }
        Aplusxiy[i+1] = (unsigned char)somme;
    }
    Aplusxiy[0]=(unsigned char)carry;
}

 void addNoOverFlowPrivateAplusxiy(unsigned int ui,__private unsigned char Aplusxiy[],__private unsigned char cteUI[],__private unsigned char PRIME[]) {

    unsigned int carry = 0;
    unsigned int somme = 0;

    multiplyNoOverFlowCte(ui,cteUI,PRIME);
    int i;
    for ( i = WORDLENGTH_T+1; i >= 1; i--)
    {
        somme = (uint)Aplusxiy[i] + (uint)cteUI[i-1]  + carry;
        carry = 0;
        if (somme >= TWOPOW_W) {
            somme -= TWOPOW_W;
            carry = 1;
        }
        Aplusxiy[i] = (unsigned char)somme;
    }
    somme = (uint)Aplusxiy[0] +  carry;
    carry = 0;
    if (somme >= TWOPOW_W) {
        somme -= TWOPOW_W;
        carry = 1;
    }
    Aplusxiy[0] = (unsigned char)somme;

}

  void multiplyNoOverFlowCte(uint n,__private unsigned char cteUI[],__private unsigned char PRIME[]) {
    int alength=WORDLENGTH_T;
    unsigned long int UV=0;
    unsigned long int U=0;
    unsigned long int V=0;
        int i;
        for( i=alength;i>=0;i--)
            cteUI[i]=0;

        for( i=alength-1;i>=0;i--) {
            U=0;
             {

                UV=(uint)cteUI[i+1] + n*((uint)PRIME[i]) + U;
                U=(UV&0xFF00)>>8;
                V=UV&0xFF;

                cteUI[i+1]=(unsigned char)V;
            }
             cteUI[i]=(unsigned char)U;
        }

    }

     void copyResultPrivTo(__global unsigned char*  outputBytes,__private unsigned char resultPrivate[] ,const size_t ID) {
        int i;
    for ( i = WORDLENGTH_T-1; i >= 0; i--) outputBytes[ID*WORDLENGTH_T+i]=    resultPrivate[i+1];

}

inline xword xaddc(xword a, xword b, xword cin, xword *cout)
{
    const uint s = (uint)a + (uint)b + (uint)cin;
    *cout = (xword)(s >> XW);
    return (xword)(s & XWMASK);
}

inline xword xsubb(xword a, xword b, xword bin, xword *bout)
{
    const int d = (int)a - (int)b - (int)bin;
    *bout = (xword)(d < 0);
    return (xword)((uint)d & XWMASK);
}

inline void xmac(xword a, xword b, xword c, xword d, xword *hi, xword *lo)
{
    const uint p = (uint)a * (uint)b + (uint)c + (uint)d;
    *lo = (xword)(p & XWMASK);
    *hi = (xword)(p >> XW);
}

inline void copyN(xword *dst, const xword *src)
{
    for (int i = 0; i < XT; i++) dst[i] = src[i];
}

inline void zeroN(xword *a)
{
    for (int i = 0; i < XT; i++) a[i] = 0;
}

inline void oneN(xword *a)
{
    zeroN(a);
    a[XT - 1] = 1;
}

inline int isZeroN(const xword *a)
{
    xword acc = 0;
    for (int i = 0; i < XT; i++) acc |= a[i];
    return acc == 0;
}

inline int cmpN(const xword *a, const xword *b)
{
    for (int i = 0; i < XT; i++) {
        if (a[i] > b[i]) return  1;
        if (a[i] < b[i]) return -1;
    }
    return 0;
}

inline xword addN(xword *r, const xword *a, const xword *b)
{
    xword carry = 0;
    for (int i = XT - 1; i >= 0; i--) r[i] = xaddc(a[i], b[i], carry, &carry);
    return carry;
}

inline xword subN(xword *r, const xword *a, const xword *b)
{
    xword borrow = 0;
    for (int i = XT - 1; i >= 0; i--) r[i] = xsubb(a[i], b[i], borrow, &borrow);
    return borrow;
}

inline xword shl1N(xword *a)
{
    xword carry = 0;
    for (int i = XT - 1; i >= 0; i--) {
        const xword nx = (xword)(a[i] >> (XW - 1));
        a[i] = (xword)((((uint)a[i] << 1) & XWMASK) | (uint)carry);
        carry = nx;
    }
    return carry;
}

inline void shr1N(xword *a)
{
    xword carry = 0;
    for (int i = 0; i < XT; i++) {
        const xword nx = (xword)(a[i] & 1u);
        a[i] = (xword)(((uint)a[i] >> 1) | ((uint)carry << (XW - 1)));
        carry = nx;
    }
}

inline xword bitAtN(const xword *a, int bit)
{
    const int w = XT - 1 - (bit / XW);
    return (xword)(((uint)a[w] >> (bit % XW)) & 1u);
}

inline int bitLenN(const xword *a)
{
    for (int i = 0; i < XT; i++) {
        if (a[i]) {
            int b = 0;
            uint v = (uint)a[i];
            while (v) { v >>= 1; b++; }
            return (XT - 1 - i) * XW + b;
        }
    }
    return 0;
}

inline void mulLowN(xword *r, const xword *a, const xword *b)
{
    xword acc[XT];
    zeroN(acc);
    for (int i = XT - 1; i >= 0; i--) {
        xword U = 0, hi, lo;
        for (int j = XT - 1; j >= 0; j--) {
            const int k = i + j + 1 - XT;
            if (k < 0) continue;
            xmac(a[i], b[j], acc[k], U, &hi, &lo);
            acc[k] = lo;
            U = hi;
        }
    }
    copyN(r, acc);
}

inline void montMulPriv(xword *out, const xword *xbe, const xword *ybe,
                        const xword *nbe, uint m_prime)
{
    xword a[XT], b[XT], n[XT], t[XT + 2];

    for (int i = 0; i < XT; i++) {
        a[i] = xbe[XT - 1 - i];
        b[i] = ybe[XT - 1 - i];
        n[i] = nbe[XT - 1 - i];
    }
    for (int i = 0; i < XT + 2; i++) t[i] = 0;

    for (int i = 0; i < XT; i++) {
        xword C = 0, hi, lo, cc;
        const xword bi = b[i];
        for (int j = 0; j < XT; j++) {
            xmac(a[j], bi, t[j], C, &hi, &lo);
            t[j] = lo; C = hi;
        }
        t[XT]     = xaddc(t[XT], C, 0, &cc);
        t[XT + 1] = cc;

        const xword m = (xword)(((uint)t[0] * m_prime) & XWMASK);
        xmac(m, n[0], t[0], 0, &hi, &lo);
        C = hi;
        for (int j = 1; j < XT; j++) {
            xmac(m, n[j], t[j], C, &hi, &lo);
            t[j - 1] = lo; C = hi;
        }
        t[XT - 1] = xaddc(t[XT], C, 0, &cc);
        t[XT]     = (xword)(((uint)t[XT + 1] + (uint)cc) & XWMASK);
    }

    int ge = (t[XT] != 0);
    if (!ge) {
        ge = 1;
        for (int j = XT - 1; j >= 0; j--) {
            if (t[j] > n[j]) { ge = 1; break; }
            if (t[j] < n[j]) { ge = 0; break; }
        }
    }
    if (ge) {
        xword borrow = 0;
        for (int j = 0; j < XT; j++) t[j] = xsubb(t[j], n[j], borrow, &borrow);
    }

    for (int i = 0; i < XT; i++) out[i] = t[XT - 1 - i];
}

inline void modDoubleN(xword *r, const xword *p)
{
    const xword c = shl1N(r);
    if (c || cmpN(r, p) >= 0) { xword tmp[XT]; subN(tmp, r, p); copyN(r, tmp); }
}

inline void reduceN(xword *r, const xword *a, const xword *p)
{
    zeroN(r);
    for (int bit = XW * XT - 1; bit >= 0; bit--) {
        const xword c = shl1N(r);
        r[XT - 1] |= bitAtN(a, bit);
        if (c || cmpN(r, p) >= 0) { xword tmp[XT]; subN(tmp, r, p); copyN(r, tmp); }
    }
}

inline void divModN(xword *q, xword *r, const xword *a, const xword *b)
{
    zeroN(q);
    zeroN(r);
    if (isZeroN(b)) return;
    for (int bit = XW * XT - 1; bit >= 0; bit--) {
        const xword c = shl1N(r);
        r[XT - 1] |= bitAtN(a, bit);
        shl1N(q);
        if (c || cmpN(r, b) >= 0) {
            xword tmp[XT];
            subN(tmp, r, b);
            copyN(r, tmp);
            q[XT - 1] |= 1;
        }
    }
}

inline void toMontN(xword *out, const xword *a, const xword *p)
{
    xword acc[XT];
    reduceN(acc, a, p);
    for (int i = 0; i < XW * XT; i++) modDoubleN(acc, p);
    copyN(out, acc);
}

inline void op_compare(__global const xword *x, __global const xword *y,
                       __global xword *out, size_t g)
{
    int c = 0;
    for (int i = 0; i < XT; i++) {
        const xword xa = x[XIDX(g, i)], yb = y[XIDX(g, i)];
        if (xa > yb) { c =  1; break; }
        if (xa < yb) { c = -1; break; }
    }
    const xword fill = (c < 0) ? (xword)XWMASK : (xword)0;
    for (int i = 0; i < XT; i++) out[XIDX(g, i)] = fill;
    out[XIDX(g, XT - 1)] = (c < 0) ? (xword)XWMASK : (xword)c;
}

inline void op_reduce(__global const xword *x, __global xword *out,
                      size_t g, const xword *p)
{
    xword a[XT], r[XT];
    for (int i = 0; i < XT; i++) a[i] = x[XIDX(g, i)];
    reduceN(r, a, p);
    for (int i = 0; i < XT; i++) out[XIDX(g, i)] = r[i];
}

inline void op_modmul(__global const xword *x, __global const xword *y,
                      __global xword *out, size_t g, const xword *p, uint m_prime)
{
    xword a[XT], b[XT], am[XT], rm[XT];
    for (int i = 0; i < XT; i++) { a[i] = x[XIDX(g, i)]; b[i] = y[XIDX(g, i)]; }
    toMontN(am, a, p);
    montMulPriv(rm, am, b, p, m_prime);
    for (int i = 0; i < XT; i++) out[XIDX(g, i)] = rm[i];
}

inline void op_modmul_r2(__global const xword *x, __global const xword *y,
                         __global xword *out, size_t g, const xword *p,
                         __constant const xword *r2g, uint m_prime)
{
    xword a[XT], b[XT], r2[XT], t1[XT], r[XT];
    for (int i = 0; i < XT; i++) {
        a[i]  = x[XIDX(g, i)];
        b[i]  = y[XIDX(g, i)];
        r2[i] = r2g[i];
    }
    montMulPriv(t1, a, b, p, m_prime);
    montMulPriv(r, t1, r2, p, m_prime);
    for (int i = 0; i < XT; i++) out[XIDX(g, i)] = r[i];
}

inline void op_modexp(__global const xword *x, __global const xword *y,
                      __global xword *out, size_t g, const xword *p, uint m_prime)
{
    xword a[XT], e[XT], am[XT], rm[XT], one[XT], r[XT];
    for (int i = 0; i < XT; i++) { a[i] = x[XIDX(g, i)]; e[i] = y[XIDX(g, i)]; }

    oneN(one);
    toMontN(am, a, p);
    toMontN(rm, one, p);

    const int nb = bitLenN(e);
    for (int bit = nb - 1; bit >= 0; bit--) {
        xword t1[XT];
        montMulPriv(t1, rm, rm, p, m_prime);
        copyN(rm, t1);
        if (bitAtN(e, bit)) {
            montMulPriv(t1, rm, am, p, m_prime);
            copyN(rm, t1);
        }
    }
    montMulPriv(r, rm, one, p, m_prime);
    for (int i = 0; i < XT; i++) out[XIDX(g, i)] = r[i];
}

inline void op_exp(__global const xword *x, __global const xword *y,
                   __global xword *out, size_t g)
{
    xword a[XT], e[XT], r[XT];
    for (int i = 0; i < XT; i++) { a[i] = x[XIDX(g, i)]; e[i] = y[XIDX(g, i)]; }
    oneN(r);
    const int nb = bitLenN(e);
    for (int bit = nb - 1; bit >= 0; bit--) {
        xword t1[XT];
        mulLowN(t1, r, r);
        copyN(r, t1);
        if (bitAtN(e, bit)) {
            mulLowN(t1, r, a);
            copyN(r, t1);
        }
    }
    for (int i = 0; i < XT; i++) out[XIDX(g, i)] = r[i];
}

inline void op_divide(__global const xword *x, __global const xword *y,
                      __global xword *out, size_t g)
{
    xword a[XT], b[XT], q[XT], r[XT];
    for (int i = 0; i < XT; i++) { a[i] = x[XIDX(g, i)]; b[i] = y[XIDX(g, i)]; }
    divModN(q, r, a, b);
    for (int i = 0; i < XT; i++) {
        out[XIDX2(g, i)]      = q[i];
        out[XIDX2(g, XT + i)] = r[i];
    }
}

inline void op_isqrt(__global const xword *x, __global xword *out, size_t g)
{
    xword a[XT], xc[XT], yc[XT], q[XT], rr[XT];
    for (int i = 0; i < XT; i++) a[i] = x[XIDX(g, i)];

    if (isZeroN(a)) {
        for (int i = 0; i < XT; i++) out[XIDX(g, i)] = 0;
        return;
    }

    const int nb = bitLenN(a);
    const int h = (nb + 1) >> 1;
    zeroN(xc);
    if (h >= XW * XT) copyN(xc, a);
    else xc[XT - 1 - (h / XW)] = (xword)(1u << (h % XW));

    for (int it = 0; it < 2 * XW * XT; it++) {
        divModN(q, rr, a, xc);
        addN(yc, xc, q);
        shr1N(yc);
        if (cmpN(yc, xc) >= 0) break;
        copyN(xc, yc);
    }
    for (int i = 0; i < XT; i++) out[XIDX(g, i)] = xc[i];
}

__kernel void mpaKernel(__global unsigned char* input1, __global unsigned char* input2, __global unsigned char* outputBytes,__constant  int* OPERATOR_WORDSIZE_BITSLENGHT_MPRIME, __constant unsigned char* globalPRIME)
{

    const uint m_prime=(uint)OPERATOR_WORDSIZE_BITSLENGHT_MPRIME[3];

    __private  unsigned char PRIME[WORDLENGTH_T];
    int i;
    for( i=0;i<WORDLENGTH_T;i++) PRIME[i]=globalPRIME[i];
    const size_t ID=get_global_id(0);
    switch(OPERATOR_WORDSIZE_BITSLENGHT_MPRIME[0]){
        case ADD : add(input1,input2,outputBytes,ID);
               break;
        case SUBTRACT : subtractPositive(input1,input2,outputBytes,ID);
                break;
        case ADDMOD : addMod(input1,input2,outputBytes,ID,PRIME);
                break;
        case SUBTRACTMOD : subtractMod(input1,input2,outputBytes,ID,PRIME);
                break;
        case MULTIPLYOPERANDSCANNING  : multiplyOperandScanning(input1,input2,outputBytes,ID);
                break;
        case MULTIPLYPRODUCTSCANNING : multiplyProductScanning(input1,input2,outputBytes,ID);
                break;
        case MONTGOMERYMULTIPLICATION :
         montgomeryMultiplication(input1,input2,outputBytes,ID,PRIME,m_prime);
               break;
        case COMPARE : op_compare(input1,input2,outputBytes,ID);
                break;
        case REDUCE : op_reduce(input1,outputBytes,ID,PRIME);
                break;
        case MODMUL : op_modmul(input1,input2,outputBytes,ID,PRIME,m_prime);
                break;
        case MODMUL_R2 : op_modmul_r2(input1,input2,outputBytes,ID,PRIME,
                                      globalPRIME+WORDLENGTH_T,m_prime);
                break;
        case MODEXP : op_modexp(input1,input2,outputBytes,ID,PRIME,m_prime);
                break;
        case EXPONENTIATION : op_exp(input1,input2,outputBytes,ID);
                break;
        case DIVIDE : op_divide(input1,input2,outputBytes,ID);
                break;
        case ISQRT : op_isqrt(input1,outputBytes,ID);
                break;
        default :
            for( i=0;i<WORDLENGTH_T;i++) outputBytes[ID*WORDLENGTH_T+i]=(unsigned char)0xFF;
        break;

    }

}
