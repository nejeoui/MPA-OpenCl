#define ADD 1
#define SUBTRACT 2
#define ADDMOD 3
#define SUBTRACTMOD 4
#define MULTIPLYOPRANDSCANNING 5
#define MULTIPLYPRODUCTSCANNING 6
#define MONTGOMERYMULTIPLICATION 7
#define TWOPOW_W 0x100000000
#include <mpaKernel_32bits.h>
void addPrime(__global uint*  outputBytes, const size_t ID, __private uint PRIME[]){
    uint carry = 0;
    int i=0;
    for ( i = WORDLENGTH_T-1; i >= 0; i--)
    {
        const size_t index=ID*WORDLENGTH_T+i;
        const ulong somme = (ulong)outputBytes[index] + (ulong)PRIME[i] + (ulong)carry;
        outputBytes[index] = (uint)somme;
        carry = (uint)(somme >> 32);
    }
}

void subtractPrime(__global uint*  outputBytes, const size_t ID,__private uint PRIME[]){
    uint borrow = 0;
    int i=0;
    for ( i = WORDLENGTH_T-1; i >= 0; i--) {
        const size_t index=ID*WORDLENGTH_T+i;
        const ulong diff = (ulong)outputBytes[index] - (ulong)PRIME[i] - (ulong)borrow;
        outputBytes[index] = (uint)diff;
        borrow = (uint)((diff >> 32) & 1u);
    }

}

char compareWithPrime(__global uint*  outputBytes, const size_t ID, __private uint PRIME[]){

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

void add(__global uint* input1, __global uint* input2, __global uint* outputBytes, const size_t ID)
{
     uint carry = 0;
     int i=0;
    for (i=WORDLENGTH_T-1; i >= 0 ; i--)
    {
        const size_t index=ID*WORDLENGTH_T+i;
        const ulong somme = (ulong)carry + (ulong)input1[index] + (ulong)input2[index];
        outputBytes[index] = (uint)somme;
        carry = (uint)(somme >> 32);
    }
}

void multiplyOperandScanning(__global uint* input1, __global uint* input2, __global uint* outputBytes, const size_t ID)
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
            UV= (ulong)outputBytes[ID*2*WORDLENGTH_T+i+j+1] + (ulong)input1[indexI]* ( (ulong)input2[indexJ])+U;
            U=(UV&0xFFFFFFFF00000000)>>32;
            V=UV&0xFFFFFFFF;

            outputBytes[ID*2*WORDLENGTH_T+ i+j+1]=(uint)V;
        }
        outputBytes[ID*2*WORDLENGTH_T+ i]=(uint)U;
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
void multiplyProductScanning(__global uint* input1, __global uint* input2, __global uint* outputBytes, const size_t ID)
{

    ulong carry = 0;
    int k;
   for( k=2*WORDLENGTH_T-2;k>=0;k--) {
        ulong accLo = carry;
        uint  accHi = 0;
        int i;
        for( i=MAX(0,k-WORDLENGTH_T+1);i<=MIN(k,WORDLENGTH_T-1);i++) {
            const ulong prod =
                ((ulong)input1[ID*WORDLENGTH_T+i]) * ((ulong)input2[ID*WORDLENGTH_T+k-i]);
            accLo += prod;
            if (accLo < prod) accHi++;
        }
        outputBytes[ID*2*WORDLENGTH_T+ k+1] = (uint)(accLo & 0xFFFFFFFF);
        carry = (accLo >> 32) | (((ulong)accHi) << 32);
    }
    outputBytes[ID*2*WORDLENGTH_T]=(uint)carry;
}

void subtractPositive(__global uint* input1, __global uint* input2, __global uint* outputBytes, const size_t ID)
{

    uint borrow = 0;
    int i=0;
    for (i = WORDLENGTH_T-1; i >= 0; i--) {
         const size_t index=ID*WORDLENGTH_T+i;
        const ulong diff = (ulong)input1[index] - (ulong)input2[index] - (ulong)borrow;
        outputBytes[index] = (uint)diff;
        borrow = (uint)((diff >> 32) & 1u);
    }
}

 void addMod(__global uint* input1, __global uint* input2, __global uint* outputBytes, const size_t ID, __private uint PRIME[])
{

    uint carry = 0;
    int i=0;

    for ( i = WORDLENGTH_T-1; i >= 0; i--)
    {
       const size_t index=ID*WORDLENGTH_T+i;
        const ulong somme = (ulong)input1[index] + (ulong)input2[index] + (ulong)carry;
        outputBytes[index] = (uint)somme;
        carry = (uint)(somme >> 32);
    }

    if (carry == 1) {

        subtractPrime(outputBytes,ID,PRIME);
    }

    else if(compareWithPrime( outputBytes,ID,PRIME)>=0) {

        subtractPrime(outputBytes,ID,PRIME);
    }

}

 void subtractMod(__global uint* input1, __global uint* input2, __global uint* outputBytes, const size_t ID,  __private uint PRIME[])
{
   uint borrow = 0;
    int i=0;
    for ( i = WORDLENGTH_T-1; i >= 0; i--)
    {
        const size_t  index=ID*WORDLENGTH_T+i;
        const ulong diff = (ulong)input1[index] - (ulong)input2[index] - (ulong)borrow;
        outputBytes[index] = (uint)diff;
        borrow = (uint)((diff >> 32) & 1u);
    }
    if (borrow == 1) {
        addPrime(outputBytes,ID,PRIME);
    }

}

  void montgomeryMultiplication(__global uint*  x,__global uint* y,__global uint* result,const size_t ID,__private uint PRIME[],const uint m_prime) {
    __private   uint resultPrivate[WORDLENGTH_T+1];
    __private   uint xiy[WORDLENGTH_T+1];
    __private   uint Aplusxiy[WORDLENGTH_T+2];

    __private   uint cteUI[WORDLENGTH_T+1];
    int i;
    const uint Yend=y[ID*WORDLENGTH_T+WORDLENGTH_T-1];
     for( i=WORDLENGTH_T;i>=0;i--) {
resultPrivate[i]=0;
            xiy[i]=0;
       Aplusxiy[i]=0;

 }

 Aplusxiy[WORDLENGTH_T+1]=0;
    for( i=WORDLENGTH_T-1;i>=0;i--) {
        size_t xindex=i+WORDLENGTH_T*ID;
        uint ui=(uint)((((ulong)resultPrivate[WORDLENGTH_T]+ ((ulong)x[xindex])*Yend)*m_prime)&0xFFFFFFFF);
        multiplyNoOverFlow1xWORDLENGTH(x[xindex],ID,y,xiy);
        addNoOverFlowPrivate_XIY(resultPrivate,xiy,Aplusxiy);
        addNoOverFlowPrivateAplusxiy(ui,Aplusxiy,cteUI,PRIME);
        rightShiftFormby1InResultPriv(Aplusxiy,resultPrivate);

    }
    if(compareResultPrivPrime(resultPrivate,PRIME)>=0) subtractPositiveResultPrivate(resultPrivate,PRIME);
    copyResultPrivTo(result,resultPrivate,ID);
}

 void rightShiftFormby1InResultPriv(__private  uint Aplusxiy[],__private  uint resultPrivate[]) {
    int i;
    for( i=0;i<WORDLENGTH_T+1;i++) resultPrivate[i]=Aplusxiy[i];
}
void subtractPositiveResultPrivate(__private uint resultPrivate[],__private uint PRIME[]){
    uint borrow = 0;
    int i;
    for ( i = WORDLENGTH_T-1; i >= 0; i--) {

        const ulong diff = (ulong)resultPrivate[i+1] - (ulong)PRIME[i] - (ulong)borrow;
        resultPrivate[i+1] = (uint)diff;
        borrow = (uint)((diff >> 32) & 1u);
    }
    resultPrivate[0]=(uint)(resultPrivate[0]-borrow);
}
int compareResultPrivPrime(__private uint resultPrivate[],__private uint PRIME[]){
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

 void multiplyNoOverFlow1xWORDLENGTH(uint n,const size_t ID,__global uint* y,__private uint xiy[]) {

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
                UV=(ulong)xiy[i+1]+ ((ulong)y[WORDLENGTH_T*ID+i])*n + U;
                 U=(UV&0xFFFFFFFF00000000)>>32;
                 V=UV&0xFFFFFFFF;

                xiy[i+1]=(uint)V;
            }
            xiy[i]=(uint)U;
        }

}

 void addNoOverFlowPrivate_XIY(__private uint resultPrivate[],__private uint xiy[],__private uint Aplusxiy[]) {
   uint carry = 0;

    int i;
    for ( i = WORDLENGTH_T; i >= 0; i--)
    {
        const ulong somme = (ulong)resultPrivate[i] + (ulong)xiy[i]  + (ulong)carry;
        Aplusxiy[i+1] = (uint)somme;
        carry = (uint)(somme >> 32);
    }
    Aplusxiy[0]=carry;
}

 void addNoOverFlowPrivateAplusxiy(uint ui,__private uint Aplusxiy[],__private uint cteUI[],__private uint PRIME[]) {

    uint carry = 0;

    multiplyNoOverFlowCte(ui,cteUI,PRIME);
    int i;
    for ( i = WORDLENGTH_T+1; i >= 1; i--)
    {
        const ulong somme = (ulong)Aplusxiy[i] + (ulong)cteUI[i-1]  + (ulong)carry;
        Aplusxiy[i] = (uint)somme;
        carry = (uint)(somme >> 32);
    }
    Aplusxiy[0] = (uint)((ulong)Aplusxiy[0] + (ulong)carry);

}

  void multiplyNoOverFlowCte(uint n,__private uint cteUI[],__private uint PRIME[]) {
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

                UV=(ulong)cteUI[i+1] + n*((ulong)PRIME[i]) + U;
                U=(UV&0xFFFFFFFF00000000)>>32;
                V=UV&0xFFFFFFFF;

                cteUI[i+1]=(uint)V;
            }
             cteUI[i]=(uint)U;
        }

    }

     void copyResultPrivTo(__global uint*  outputBytes,__private uint resultPrivate[] ,const size_t ID) {
        int i;
    for ( i = WORDLENGTH_T-1; i >= 0; i--) outputBytes[ID*WORDLENGTH_T+i]=    resultPrivate[i+1];

}

__kernel void mpaKernel(__global uint* input1, __global uint* input2, __global uint* outputBytes,__constant  int* OPERATOR_WORDSIZE_BITSLENGHT_MPRIME, __constant uint* globalPRIME)
{

    const uint m_prime=(uint)OPERATOR_WORDSIZE_BITSLENGHT_MPRIME[3];

    __private  uint PRIME[WORDLENGTH_T];
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
        case MULTIPLYOPRANDSCANNING : multiplyOperandScanning(input1,input2,outputBytes,ID);
                break;
        case MULTIPLYPRODUCTSCANNING : multiplyProductScanning(input1,input2,outputBytes,ID);
                break;
        case MONTGOMERYMULTIPLICATION :
         montgomeryMultiplication(input1,input2,outputBytes,ID,PRIME,m_prime);
        break;
        default :
            for( i=0;i<WORDLENGTH_T;i++) outputBytes[ID*WORDLENGTH_T+i]=0xFFFFFFFF;
        break;

    }

}
