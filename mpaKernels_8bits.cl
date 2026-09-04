#define TWOPOW_W 256
#define ADD 1
#define SUBTRACT 2
#define ADDMOD 3
#define SUBTRACTMOD 4
#define MULTIPLYOPRANDSCANNING 5
#define MULTIPLYPRODUCTSCANNING 6
#define MONTGOMERYMULTIPLICATION 7
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
        case MULTIPLYOPRANDSCANNING  : multiplyOperandScanning(input1,input2,outputBytes,ID);
                break;
        case MULTIPLYPRODUCTSCANNING : multiplyProductScanning(input1,input2,outputBytes,ID);
                break;
        case MONTGOMERYMULTIPLICATION :
         montgomeryMultiplication(input1,input2,outputBytes,ID,PRIME,m_prime);
               break;
        default :
            for( i=0;i<WORDLENGTH_T;i++) outputBytes[ID*WORDLENGTH_T+i]=(unsigned char)0xFF;
        break;

    }

}
