#define TWOPOW_W 65536
#define ADD 1
#define SUBTRACT 2
#define ADDMOD 3
#define SUBTRACTMOD 4
#define MULTIPLYOPRANDSCANNING 5
#define MULTIPLYPRODUCTSCANNING 6
#define MONTGOMERYMULTIPLICATION 7
#include <mpaKernel_16bits.h>

// use mul_hi 
void addPrime(__global ushort*  outputBytes, const size_t ID, __private ushort PRIME[]){
    uint carry = 0;
    uint somme = 0;
    int i=0;
    for ( i = WORDLENGTH_T-1; i >= 0; i--) 
    {
        const size_t index=ID*WORDLENGTH_T+i;
        somme = carry+(uint)outputBytes[index] + (uint)PRIME[i] ;
        carry = 0;
        if (somme >= TWOPOW_W) {
            somme -= TWOPOW_W;
            carry = 1;
        }
        outputBytes[index] =(ushort) somme;
    }
}

void subtractPrime(__global ushort*  outputBytes, const size_t ID,__private ushort PRIME[]){
    int borrow = 0;
    int diff = 0;
    int i=0;
    for ( i = WORDLENGTH_T-1; i >= 0; i--) {
        const size_t index=ID*WORDLENGTH_T+i;
        diff = (uint)outputBytes[index] - (uint)PRIME[i] - borrow;
        borrow = 0;
        if (diff < 0) {
            diff += TWOPOW_W;
            borrow = 1;
        }
        outputBytes[index] = (ushort)diff;
    }

}

char compareWithPrime(__global ushort*  outputBytes, const size_t ID, __private ushort PRIME[]){
    
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


void add(__global ushort* input1, __global ushort* input2, __global ushort* outputBytes, const size_t ID)
{
    uint carry = 0;
    int somme = 0;
    int i=0;
    for (i=WORDLENGTH_T-1; i >= 0 ; i--)
    {
        const size_t index=ID*WORDLENGTH_T+i;
        somme = carry + (uint)input1[index] + (uint)input2[index];
        carry = 0;
        if (somme >= TWOPOW_W) {
            somme -= TWOPOW_W;
            carry = 1;
        }
        outputBytes[index] =(ushort) somme;
        
    }
}

void multiplyOperandScanning(__global ushort* input1, __global ushort* input2, __global ushort* outputBytes, const size_t ID)
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
            UV= (uint)outputBytes[ID*2*WORDLENGTH_T+i+j+1]+ ((uint)input1[indexI])*input2[indexJ]+U;
            U=(UV&0xFFFF0000)>>16;
            V=UV&0xFFFF;
            
            outputBytes[ID*2*WORDLENGTH_T+ i+j+1]=(ushort)V;
        }
        outputBytes[ID*2*WORDLENGTH_T+ i]=(ushort)U;
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
void multiplyProductScanning(__global ushort* input1, __global ushort* input2, __global ushort* outputBytes, const size_t ID)
{
    unsigned long int UV=0;
    unsigned long int U=0;
    unsigned long int V=0;
    int k;
   for( k=2*WORDLENGTH_T-2;k>=0;k--) {
        UV=0;
        int i;
        for( i=MAX(0,k-WORDLENGTH_T+1);i<=MIN(k,WORDLENGTH_T-1);i++) 
            UV+=input1[ID*WORDLENGTH_T+i]*input2[ID*WORDLENGTH_T+k-i];
        UV=UV+U;
        U=(UV&0xFFFF0000)>>16;
        V=UV&0xFFFF;
        outputBytes[ID*2*WORDLENGTH_T+ k+1]=(ushort)V;
        
    }
    outputBytes[ID*2*WORDLENGTH_T]=(ushort)U;
}
 
void subtractPositive(__global ushort* input1, __global ushort* input2, __global ushort* outputBytes, const size_t ID)
{
    
    int borrow = 0;
    int diff = 0;
    int i=0;
    for (i = WORDLENGTH_T-1; i >= 0; i--) {
         const size_t index=ID*WORDLENGTH_T+i;
        diff = (uint)input1[index] - (uint)input2[index] - borrow;
        borrow = 0;
        if (diff < 0) {
            diff += TWOPOW_W;
            borrow = 1;
        }
        outputBytes[index] =(ushort) diff;
    }
    
   
}

 void addMod(__global ushort* input1, __global ushort* input2, __global ushort* outputBytes, const size_t ID, __private ushort PRIME[])
{
   
   
    int carry = 0;
    int somme = 0;
    int i=0;

    for ( i = WORDLENGTH_T-1; i >= 0; i--) 
    {
       const size_t index=ID*WORDLENGTH_T+i;
        somme = (uint)input1[index] + (uint)input2[index] + carry;
        carry = 0;
        if (somme >= TWOPOW_W) {
            somme -= TWOPOW_W;
            carry = 1;
        }
        outputBytes[index] =(ushort) somme;
    }
    
    if (carry == 1) {
        
        subtractPrime(outputBytes,ID,PRIME);
    }
    
    else if(compareWithPrime( outputBytes,ID,PRIME)==1) {
        
        subtractPrime(outputBytes,ID,PRIME);
    }
    
    
}



 void subtractMod(__global ushort* input1, __global ushort* input2, __global ushort* outputBytes, const size_t ID,  __private ushort PRIME[])
{
   int borrow = 0;
    int diff = 0;
    int i=0;
    for ( i = WORDLENGTH_T-1; i >= 0; i--) 
    {
        const size_t  index=ID*WORDLENGTH_T+i;
        diff = (uint)input1[index] - (uint)input2[index] - borrow;
        borrow = 0;
        if (diff < 0) {
            diff += TWOPOW_W;
            borrow = 1;
        }
        
        outputBytes[index] =(ushort) diff;
    }
    if (borrow == 1) {
        addPrime(outputBytes,ID,PRIME);
    }
    
    
} 

  void montgomeryMultiplication(__global ushort*  x,__global ushort* y,__global ushort* result,const size_t ID,__private ushort PRIME[],const unsigned int m_prime) {
    __private   ushort resultPrivate[WORDLENGTH_T+1];
    __private   ushort xiy[WORDLENGTH_T+1];
    __private   ushort Aplusxiy[WORDLENGTH_T+2];
    

    __private   ushort cteUI[WORDLENGTH_T+1];
    int i;
    const ushort Yend=y[ID*WORDLENGTH_T+WORDLENGTH_T-1];
     for( i=WORDLENGTH_T;i>=0;i--) {
resultPrivate[i]=0;
            xiy[i]=0;
       Aplusxiy[i]=0;
          
 }
 
 Aplusxiy[WORDLENGTH_T+1]=0;
    for( i=WORDLENGTH_T-1;i>=0;i--) {
        size_t xindex=i+WORDLENGTH_T*ID;
        unsigned short ui=(ushort)(((uint)resultPrivate[WORDLENGTH_T]+ ((uint)x[xindex])*Yend)*m_prime);
        multiplyNoOverFlow1xWORDLENGTH(x[xindex],ID,y,xiy);
        addNoOverFlowPrivate_XIY(resultPrivate,xiy,Aplusxiy);
        addNoOverFlowPrivateAplusxiy(ui,Aplusxiy,cteUI,PRIME);
        rightShiftFormby1InResultPriv(Aplusxiy,resultPrivate);
        
    }
    if(compareResultPrivPrime(resultPrivate,PRIME)==1) subtractPositiveResultPrivate(resultPrivate,PRIME);
    copyResultPrivTo(result,resultPrivate,ID);
}

 void rightShiftFormby1InResultPriv(__private  ushort Aplusxiy[],__private  ushort resultPrivate[]) {
    int i;
    for( i=0;i<WORDLENGTH_T+1;i++) resultPrivate[i]=Aplusxiy[i];
}
void subtractPositiveResultPrivate(__private ushort resultPrivate[],__private ushort PRIME[]){
    int borrow = 0;
    int diff = 0;
    int i;
    for ( i = WORDLENGTH_T-1; i >= 0; i--) {
        
        diff = (uint)resultPrivate[i+1] - (uint)PRIME[i] - borrow;
        borrow = 0;
        if (diff < 0) {
            diff += TWOPOW_W;
            borrow = 1;
        }
        resultPrivate[i+1] = (ushort)diff;
    }
    resultPrivate[0]=(ushort)(resultPrivate[0]-borrow);
}
int compareResultPrivPrime(__private ushort resultPrivate[],__private ushort PRIME[]){
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

 void multiplyNoOverFlow1xWORDLENGTH(ushort n,const size_t ID,__global ushort* y,__private ushort xiy[]) {
    
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
                U=(UV&0x00000000FFFF0000)>>16;
                V=(ushort)UV;
                
                xiy[i+1]=(ushort)V;
            }
            xiy[i]=(ushort)U;
        }
    
}

 void addNoOverFlowPrivate_XIY(__private ushort resultPrivate[],__private ushort xiy[],__private ushort Aplusxiy[]) {
   unsigned int carry = 0;
   unsigned  int somme = 0;

    int i;
    for ( i = WORDLENGTH_T; i >= 0; i--) 
    {
        somme = (uint)resultPrivate[i] + (uint)xiy[i]  + carry;
        carry = 0;
        if (somme >= TWOPOW_W) {
            somme %= TWOPOW_W;
            carry = 1;
        }
        Aplusxiy[i+1] = (ushort)somme;
    }
    Aplusxiy[0]=(ushort)carry;
}

 void addNoOverFlowPrivateAplusxiy(unsigned int ui,__private ushort Aplusxiy[],__private ushort cteUI[],__private ushort PRIME[]) {
    
    unsigned int carry = 0;
    unsigned int somme = 0;
    
    multiplyNoOverFlowCte(ui,cteUI,PRIME);
    int i;
    for ( i = WORDLENGTH_T+1; i >= 1; i--) 
    {
        somme = (uint)Aplusxiy[i] + (uint)cteUI[i-1]  + carry;
        carry = 0;
        if (somme >= TWOPOW_W) {
            somme %= TWOPOW_W;
            carry = 1;
        }
        Aplusxiy[i] = (ushort)somme;
    }
    somme = Aplusxiy[0] +  carry;
    carry = 0;
    if (somme >= TWOPOW_W) {
        somme %= TWOPOW_W;
        carry = 1;
    }
    Aplusxiy[0] = (ushort)somme;
    
}

  void multiplyNoOverFlowCte(int n,__private ushort cteUI[],__private ushort PRIME[]) {
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
                U=(UV&0x00000000FFFF0000)>>16;
                V=(ushort)UV;
                

                cteUI[i+1]=(ushort)V;
            }
             cteUI[i]=(ushort)U;
        }
        
    }

     void copyResultPrivTo(__global ushort*  outputBytes,__private ushort resultPrivate[] ,const size_t ID) {
        int i;
    for ( i = WORDLENGTH_T-1; i >= 0; i--) outputBytes[ID*WORDLENGTH_T+i]=    resultPrivate[i+1];
    
}

__kernel void mpaKernel(__global ushort* input1, __global ushort* input2, __global ushort* outputBytes,__constant  int* OPERATOR_WORDSIZE_BITSLENGHT_MPRIME, __constant ushort* globalPRIME)
{
    
    const int m_prime=OPERATOR_WORDSIZE_BITSLENGHT_MPRIME[3];
    
    __private  ushort PRIME[16]={0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xfffe,0xffff,0xfc2f,};
   // for(int i=0;i<WORDLENGTH_T;i++) PRIME[i]=globalPRIME[i];
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

