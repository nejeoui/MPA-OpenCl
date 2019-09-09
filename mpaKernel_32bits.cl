
#define ADD 1
#define SUBTRACT 2
#define ADDMOD 3
#define SUBTRACTMOD 4
#define MULTIPLYOPRANDSCANNING 5
#define MULTIPLYPRODUCTSCANNING 6
#define MONTGOMERYMULTIPLICATION 7
#define TWOPOW_W 0x100000000 
#include <mpaKernel_32bits.h>
// use mul_hi 
void addPrime(__global uint*  outputBytes, const size_t ID, __private uint PRIME[]){
    char carry = 0;
    ulong somme = 0;
    int i=0;
    for ( i = WORDLENGTH_T-1; i >= 0; i--) 
    {
        const size_t index=ID*WORDLENGTH_T+i;
        somme = outputBytes[index] + PRIME[i] + carry;
        carry = 0;
        if (somme >= TWOPOW_W) {
            somme %= TWOPOW_W;
            carry = 1;
        }
        outputBytes[index] =(uint) somme;
    }
}

void subtractPrime(__global uint*  outputBytes, const size_t ID,__private uint PRIME[]){
    long borrow = 0;
    long diff = 0;
    int i=0;
    for ( i = WORDLENGTH_T-1; i >= 0; i--) {
        const size_t index=ID*WORDLENGTH_T+i;
        diff = (ulong)outputBytes[index] - (ulong)PRIME[i] - borrow;
        borrow = 0;
        if (diff < 0) {
            diff += TWOPOW_W;
            borrow = 1;
        }
        outputBytes[index] = (uint)diff;
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
     ulong carry = 0;
     ulong somme = 0;
     int i=0;
    for (i=WORDLENGTH_T-1; i >= 0 ; i--)
    {
        const size_t index=ID*WORDLENGTH_T+i;
        somme = carry+ (ulong)input1[index] + (ulong)input2[index] ;
        carry = 0;
        if (somme >= TWOPOW_W) {
            somme -= TWOPOW_W;
            carry = 1;
        }
        outputBytes[index] =(uint) somme;
        
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
    unsigned long int UV=0;
    unsigned long int U=0;
    unsigned long int V=0;
    int k;
   for( k=2*WORDLENGTH_T-2;k>=0;k--) {
        UV=0;
        int i;
        for( i=MAX(0,k-WORDLENGTH_T+1);i<=MIN(k,WORDLENGTH_T-1);i++) 
            UV+=((ulong)input1[ID*WORDLENGTH_T+i])* ((ulong)input2[ID*WORDLENGTH_T+k-i]);
        UV=UV+U;
        U=(UV&0xFFFFFFFF00000000)>>32;
        V=UV&0xFFFFFFFF;
        outputBytes[ID*2*WORDLENGTH_T+ k+1]=(uint)V;
        
    }
    outputBytes[ID*2*WORDLENGTH_T]=(uint)U;
}
 
void subtractPositive(__global uint* input1, __global uint* input2, __global uint* outputBytes, const size_t ID)
{
    
    long borrow = 0;
    long diff = 0;
    int i=0;
    for (i = WORDLENGTH_T-1; i >= 0; i--) {
         const size_t index=ID*WORDLENGTH_T+i;
        diff = (ulong)input1[index] - (ulong)input2[index] - borrow;
        borrow = 0;
        if (diff < 0) {
            diff += TWOPOW_W;
            borrow = 1;
        }
        outputBytes[index] =(uint) diff;
    }
    
   
}

 void addMod(__global uint* input1, __global uint* input2, __global uint* outputBytes, const size_t ID, __private uint PRIME[])
{
   
   
    ulong carry = 0;
    ulong somme = 0;
    int i=0;

    for ( i = WORDLENGTH_T-1; i >= 0; i--) 
    {
       const size_t index=ID*WORDLENGTH_T+i;
        somme = (ulong)input1[index] + (ulong)input2[index] + carry;
        carry = 0;
        if (somme >= TWOPOW_W) {
            somme -= TWOPOW_W;
            carry = 1;
        }
        outputBytes[index] =(uint) somme;
    }
    
    if (carry == 1) {
        
        subtractPrime(outputBytes,ID,PRIME);
    }
    
    else if(compareWithPrime( outputBytes,ID,PRIME)==1) {
        
        subtractPrime(outputBytes,ID,PRIME);
    }
    
    
}



 void subtractMod(__global uint* input1, __global uint* input2, __global uint* outputBytes, const size_t ID,  __private uint PRIME[])
{
   long borrow = 0;
    long diff = 0;
    int i=0;
    for ( i = WORDLENGTH_T-1; i >= 0; i--) 
    {
        const size_t  index=ID*WORDLENGTH_T+i;
        diff = (ulong)input1[index] - (ulong)input2[index] - borrow;
        borrow = 0;
        if (diff < 0) {
            diff += TWOPOW_W;
            borrow = 1;
        }
        
        outputBytes[index] =(uint) diff;
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
        ulong ui=(((ulong)resultPrivate[WORDLENGTH_T]+ ((ulong)x[xindex])*Yend)*m_prime)&0xFFFFFFFF;
        multiplyNoOverFlow1xWORDLENGTH(x[xindex],ID,y,xiy);
        addNoOverFlowPrivate_XIY(resultPrivate,xiy,Aplusxiy);
        addNoOverFlowPrivateAplusxiy(ui,Aplusxiy,cteUI,PRIME);
        rightShiftFormby1InResultPriv(Aplusxiy,resultPrivate);
        
    }
    if(compareResultPrivPrime(resultPrivate,PRIME)==1) subtractPositiveResultPrivate(resultPrivate,PRIME);
    copyResultPrivTo(result,resultPrivate,ID);
}

 void rightShiftFormby1InResultPriv(__private  uint Aplusxiy[],__private  uint resultPrivate[]) {
    int i;
    for( i=0;i<WORDLENGTH_T+1;i++) resultPrivate[i]=Aplusxiy[i];
}
void subtractPositiveResultPrivate(__private uint resultPrivate[],__private uint PRIME[]){
    long borrow = 0;
    long diff = 0;
    int i;
    for ( i = WORDLENGTH_T-1; i >= 0; i--) {
        
        diff = (ulong)resultPrivate[i+1] - (ulong)PRIME[i] - borrow;
        borrow = 0;
        if (diff < 0) {
            diff += TWOPOW_W;
            borrow = 1;
        }
        resultPrivate[i+1] = (uint)diff;
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
   ulong carry = 0;
   ulong somme = 0;

    int i;
    for ( i = WORDLENGTH_T; i >= 0; i--) 
    {
        somme = (ulong)resultPrivate[i] + (ulong)xiy[i]  + carry;
        carry = 0;
        if (somme >= TWOPOW_W) {
            somme -= TWOPOW_W;
            carry = 1;
        }
        Aplusxiy[i+1] = (uint)somme;
    }
    Aplusxiy[0]=(uint)carry;
}

 void addNoOverFlowPrivateAplusxiy(ulong ui,__private uint Aplusxiy[],__private uint cteUI[],__private uint PRIME[]) {
    
    ulong carry = 0;
    ulong somme = 0;
    
    multiplyNoOverFlowCte(ui,cteUI,PRIME);
    int i;
    for ( i = WORDLENGTH_T+1; i >= 1; i--) 
    {
        somme = (ulong)Aplusxiy[i] + (ulong)cteUI[i-1]  + carry;
        carry = 0;
        if (somme >= TWOPOW_W) {
            somme -= TWOPOW_W;
            carry = 1;
        }
        Aplusxiy[i] = (uint)somme;
    }
    somme = Aplusxiy[0] +  carry;
    carry = 0;
    if (somme >= TWOPOW_W) {
        somme -= TWOPOW_W;
        carry = 1;
    }
    Aplusxiy[0] = (uint)somme;
    
}

  void multiplyNoOverFlowCte(int n,__private uint cteUI[],__private uint PRIME[]) {
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
    
    const int m_prime=OPERATOR_WORDSIZE_BITSLENGHT_MPRIME[3];
    
    __private  uint PRIME[WORDLENGTH_T]
     ={0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xfffffffe,0xfffffc2f};
  //  int i;
 //   for( i=0;i<WORDLENGTH_T;i++) PRIME[i]=globalPRIME[i];
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
        

    }
    

}  
