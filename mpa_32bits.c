// add Unit test using cmocka.org
/*
https://community.amd.com/thread/145960
Okay, I see what the problem is now.
I had 'uint64_t' defined as 'unsigned long long' (which is standard in 32-bit gcc and also works in 64-bit gcc).
But in the OpenCL spec, the 64-bit type is 'unsigned long', and 'unsigned long long' is the 128-bit type.
If I define uint64_t as 'ulong' as per the OpenCL spec, mul_hi works correctly.
The reason why didn't blow up earlier is that, apparently, Stream does NOT really treat 'unsigned long long' as a 128-bit type. (So it does not get off scot-free, there's still a bug there). In particular, sizeof(unsigned long long) is 8,
 and all my code except for the mul_hi instruction works as if it were uint64_t.

 64-bit Atomics
 
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include <openssl/rand.h>
#include <time.h>
#include <limits.h>
#include <errno.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif
#define COMPARE 0
#define ADD 1
#define SUBTRACT 2
#define ADDMOD 3
#define SUBTRACTMOD 4
#define MULTIPLYOPRANDSCANNING 5
#define MULTIPLYPRODUCTSCANNING 6
#define MONTGOMERYMULTIPLICATION 7
#define ARITHMETICS 8
#define MAX_SOURCE_SIZE (0x100000)
// Terminal Colors
#define  END_COLOR   "\x1b[0m"
#define  BLUE_TERMINAL    "\x1b[34m"
#define  RED_TERMINAL     "\x1b[31m"
#define  GREEN_TERMINAL   "\x1b[32m"


int compareArray(unsigned int*  a,unsigned int*  b,const int SIZE,const int WORDLINGTH) {
    int diff=WORDLINGTH-SIZE;
for (int i = 0; i < SIZE; i++){
if(a[i+diff]>b[i]) {
//printf(">  %d %d %d\n",i,a[i],b[i]);
    return 1;}
if(a[i+diff]<b[i]) {
//printf("<  %d %d %d\n",i,a[i],b[i]);
    return -1;
}
}
return 0;
}
char ansi[]= "\x1b[32m";
void printArray(unsigned int*  bytes,const int SIZE,const size_t ID) {
  printf("[\n");
   for (int i = 0; i < SIZE-1; i++){
     printf("%u,",bytes[ID*SIZE+i]);
     }
     printf("%u]\n",bytes[ID*SIZE+SIZE-1]);
     
    }


void printDeviceInfo(cl_device_id device)
{   
    char queryBuffer[1024];
    int queryInt;
    cl_int clError;
    clError = clGetDeviceInfo(device, CL_DEVICE_NAME,
                              sizeof(queryBuffer),
                              &queryBuffer, NULL);
    printf("CL_DEVICE_NAME: %s\n", queryBuffer);
    queryBuffer[0] = '\0';
    clError = clGetDeviceInfo(device, CL_DEVICE_VENDOR,
                              sizeof(queryBuffer), &queryBuffer,
                              NULL);
    printf("CL_DEVICE_VENDOR: %s\n", queryBuffer);
    queryBuffer[0] = '\0';
    clError = clGetDeviceInfo(device, CL_DRIVER_VERSION,
                              sizeof(queryBuffer), &queryBuffer,
                              NULL);
    printf("CL_DRIVER_VERSION: %s\n", queryBuffer);
    queryBuffer[0] = '\0';
    clError = clGetDeviceInfo(device, CL_DEVICE_VERSION,
                              sizeof(queryBuffer), &queryBuffer,
                              NULL);
    printf("CL_DEVICE_VERSION: %s\n", queryBuffer);
    queryBuffer[0] = '\0';
    clError = clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS,
                              sizeof(int), &queryInt, NULL);
    printf("CL_DEVICE_MAX_COMPUTE_UNITS: %d\n", queryInt);
}
const char *decode(int OPERATOR){
switch(OPERATOR){
        case ADD: return "ADD"; 
        case ADDMOD: return "ADDMOD";
        case SUBTRACTMOD: return "SUBTRACTMOD";
        case SUBTRACT: return "SUBTRACT";
        case MULTIPLYPRODUCTSCANNING: return "MULTIPLYPRODUCTSCANNING";
        case MULTIPLYOPRANDSCANNING: return "MULTIPLYOPRANDSCANNING";
        case MONTGOMERYMULTIPLICATION: return "MONTGOMERYMULTIPLICATION";
       
    }
    return "INDEFINED OPERATOR";
}

const char *getErrorString(cl_int error)
{
    switch(error){
        
        case 0: return "CL_SUCCESS";
        case -1: return "CL_DEVICE_NOT_FOUND";
        case -2: return "CL_DEVICE_NOT_AVAILABLE";
        case -3: return "CL_COMPILER_NOT_AVAILABLE";
        case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
        case -5: return "CL_OUT_OF_RESOURCES";
        case -6: return "CL_OUT_OF_HOST_MEMORY";
        case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
        case -8: return "CL_MEM_COPY_OVERLAP";
        case -9: return "CL_IMAGE_FORMAT_MISMATCH";
        case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
        case -11: return "CL_BUILD_PROGRAM_FAILURE";
        case -12: return "CL_MAP_FAILURE";
        case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
        case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
        case -15: return "CL_COMPILE_PROGRAM_FAILURE";
        case -16: return "CL_LINKER_NOT_AVAILABLE";
        case -17: return "CL_LINK_PROGRAM_FAILURE";
        case -18: return "CL_DEVICE_PARTITION_FAILED";
        case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";
        case -30: return "CL_INVALID_VALUE";
        case -31: return "CL_INVALID_DEVICE_TYPE";
        case -32: return "CL_INVALID_PLATFORM";
        case -33: return "CL_INVALID_DEVICE";
        case -34: return "CL_INVALID_CONTEXT";
        case -35: return "CL_INVALID_QUEUE_PROPERTIES";
        case -36: return "CL_INVALID_COMMAND_QUEUE";
        case -37: return "CL_INVALID_HOST_PTR";
        case -38: return "CL_INVALID_MEM_OBJECT";
        case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
        case -40: return "CL_INVALID_IMAGE_SIZE";
        case -41: return "CL_INVALID_SAMPLER";
        case -42: return "CL_INVALID_BINARY";
        case -43: return "CL_INVALID_BUILD_OPTIONS";
        case -44: return "CL_INVALID_PROGRAM";
        case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
        case -46: return "CL_INVALID_KERNEL_NAME";
        case -47: return "CL_INVALID_KERNEL_DEFINITION";
        case -48: return "CL_INVALID_KERNEL";
        case -49: return "CL_INVALID_ARG_INDEX";
        case -50: return "CL_INVALID_ARG_VALUE";
        case -51: return "CL_INVALID_ARG_SIZE";
        case -52: return "CL_INVALID_KERNEL_ARGS";
        case -53: return "CL_INVALID_WORK_DIMENSION";
        case -54: return "CL_INVALID_WORK_GROUP_SIZE";
        case -55: return "CL_INVALID_WORK_ITEM_SIZE";
        case -56: return "CL_INVALID_GLOBAL_OFFSET";
        case -57: return "CL_INVALID_EVENT_WAIT_LIST";
        case -58: return "CL_INVALID_EVENT";
        case -59: return "CL_INVALID_OPERATION";
        case -60: return "CL_INVALID_GL_OBJECT";
        case -61: return "CL_INVALID_BUFFER_SIZE";
        case -62: return "CL_INVALID_MIP_LEVEL";
        case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
        case -64: return "CL_INVALID_PROPERTY";
        case -65: return "CL_INVALID_IMAGE_DESCRIPTOR";
        case -66: return "CL_INVALID_COMPILER_OPTIONS";
        case -67: return "CL_INVALID_LINKER_OPTIONS";
        case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";
        case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
        case -1001: return "CL_PLATFORM_NOT_FOUND_KHR";
        case -1002: return "CL_INVALID_D3D10_DEVICE_KHR";
        case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR";
        case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
        case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
        default: return "Unknown OpenCL error";
    }
}

int main(int argc, char **argv)
{
 int base=10;
 int DEBUG_MODE=0;
 char *endptr, *str;
unsigned long long int iterations;

           if (argc < 5) {
               fprintf(stderr, "Usage: %s NumberOfIteration OPERATOR BITLENGTH WORDSIZE [DEBUG_MODE] \n", argv[0]);
               exit(EXIT_FAILURE);
           }

           str = argv[1];

           errno = 0;    /* To distinguish success/failure after call */
           iterations = strtoull(str, &endptr, base);

           /* Check for various possible errors */

           if ((errno == ERANGE && (iterations == LONG_MAX || iterations == LONG_MIN))
                   || (errno != 0 && iterations == 0)) {
               perror("strtoull");
               exit(EXIT_FAILURE);
           }

           if (endptr == str) {
               fprintf(stderr, "Please specify these arguments NumberOfIteration OPERATOR BITLENGTH WORDSIZE  !\n");
               exit(EXIT_FAILURE);
           }

           DEBUG_MODE = (argc > 5) ? atoi(argv[5]) : 0;
    
    int OPERATOR = atoi(argv[2]);

    int WORDSIZE = atoi(argv[4]);

     int BITSLENGTH = atoi(argv[3]);

     int WORDLENGTH = BITSLENGTH/WORDSIZE;  
   //  printf("WORDLENGTH = %d\n",WORDLENGTH);
       int MPRIME = 49;

    unsigned int* PRIME;
    PRIME= (unsigned int*)malloc(WORDLENGTH*sizeof(unsigned int));
     mpz_t bigPrime;
    switch(BITSLENGTH){
        case 256 : { 
        const char* primeStr= "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";
         
          mpz_init(bigPrime);
          size_t* count;
          count = (size_t*) malloc(sizeof(size_t));
          mpz_set_ui(bigPrime,0);
          mpz_set_str(bigPrime,primeStr, 16);

          mpz_export(PRIME, count, 1, WORDLENGTH*sizeof(unsigned int), 1, 0, bigPrime);

          //precalculated m_prime
            MPRIME = 49;
            }
             break;
        case 512 : { 
          //  brainpoolP512r1 https://www.teletrust.de/fileadmin/files/oid/ecgdsa_final.pdf
        const char* primeStr=  
          "AADD9DB8DBE9C48B3FD4E6AE33C9FC07CB308DB3B3C9D20ED6639CCA703308717D4D9B009BC66842AECDA12AE6A380E62881FF2F2D82C68528AA6056583A48F3";
          
          mpz_init(bigPrime);
          size_t* count;
          count = (size_t*) malloc(sizeof(size_t));
          mpz_set_ui(bigPrime,0);
          mpz_set_str(bigPrime,primeStr, 16);
          mpz_export(PRIME, count, 1, WORDLENGTH*sizeof(unsigned int), 1, 0, bigPrime);

          //precalculated m_prime
            MPRIME = 49;
            }
             break;
        case 1024 : { 
        const char* primeStr= 
          "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"\
          "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"\
          "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"\
          "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF97";
          mpz_init(bigPrime);
          size_t* count;
          count = (size_t*) malloc(sizeof(size_t));
          mpz_set_ui(bigPrime,0);
          mpz_set_str(bigPrime,primeStr, 16);
          mpz_export(PRIME, count, 1, WORDLENGTH*sizeof(unsigned int), 1, 0, bigPrime);

          //precalculated m_prime
            MPRIME = 49;
            }
             break;
        case 2048 :{ 
        const char* primeStr= 
          "E53DF5FC3F650D066875837012A4E7BEA863C65CB592D9C36942CF69CBC6DD4F"\
          "D804E19CCF2696C9BEBCF18742FA5FB091CBDE1782E8291009464913ECE19745"\
          "7800EA6E43B0E2A64615D182B6DE150479C58D1C7C702D47EA3031B379CA13A2"\
          "048C964E1D1E8D4CD3815D0895BF31E53271D4607E16461B77FB26100915D679"\
          "9060203EDEBFEA9495A5A8E7CED68FC9DB2D47CE7992461BA78174608AD0BBE3"\
          "F5E63EC6C960564430CBD2E6E587D08EE12F94B5B99DFFB12C6727A25E800DAC"\
          "6CD8DE77A5BBC93B36E444B070888CB5ADD991870466968A6E9A23C2EE0A1D67"\
          "1C9B601081A44AA6A58D4DC76686EF15FCE1C9AEB4033395A9B24BE1AA1929BB";
          mpz_init(bigPrime);
          size_t* count;
          count = (size_t*) malloc(sizeof(size_t));
          mpz_set_ui(bigPrime,0);
          mpz_set_str(bigPrime,primeStr, 16);
          mpz_export(PRIME, count, 1, WORDLENGTH*sizeof(unsigned int), 1, 0, bigPrime);
          //precalculated m_prime
            MPRIME = 49;
            }
             break;
    }
    
 //  printArray(PRIME,WORDLENGTH,0);
    struct timespec tstart={0,0}, tend_init={0,0} , tend_createContext={0,0}, 
    tend_loadTomemory={0,0},tend_BuildProgram={0,0}, tend_createKernel={0,0}, tend_exec={0,0}, tend_redResults={0,0}, tend_test={0,0};
    clock_gettime(CLOCK_MONOTONIC, &tstart);

    cl_platform_id* platform_id = NULL;
    platform_id=(cl_platform_id*) malloc(sizeof(cl_platform_id) * 4);
    cl_device_id device_id = NULL;
    cl_context context = NULL;
    cl_command_queue command_queue = NULL;
    cl_mem Amobj = NULL;
    cl_mem Bmobj = NULL;
    cl_mem Cmobj = NULL;
    cl_mem Omobj = NULL;
    cl_mem Pmobj = NULL;
    cl_program program = NULL;
    cl_kernel kernel = NULL;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret;
    size_t k1=1024;
    size_t K = k1*iterations;
  //  printf("K = %zu\n", K);
   const size_t global[]={K/WORDLENGTH}; // global domain size
   const size_t local[]={1};
    int i, j;
    unsigned int* A;
    unsigned int* B;
    unsigned int* C;
    int* OPERATOR_WORDSIZE_BITSLENGHT_MPRIME;
    
    A = (unsigned int*)malloc(K*sizeof(unsigned int));
    B = (unsigned int*)malloc(K*sizeof(unsigned int));
    if(OPERATOR==MULTIPLYOPRANDSCANNING||OPERATOR==MULTIPLYPRODUCTSCANNING) C = (unsigned int*)malloc(2*K*sizeof(unsigned int));
    else C = (unsigned int*)malloc(K*sizeof(unsigned int));
    OPERATOR_WORDSIZE_BITSLENGHT_MPRIME = (int*)malloc(4*sizeof(int));
    
    
    FILE *fp;
    const char fileName[] = "mpaKernel_32bits.cl";
    size_t source_size;
     char *source_str;
    
    /* Load kernel source file */
    fp = fopen(fileName, "rb");
    if (!fp) {
        fprintf(stderr, "Failed to load kernel.\n");
        exit(1);
    }
    source_str = ( char *)malloc(MAX_SOURCE_SIZE);
    source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
    fclose(fp);
    
    /* Initialize input data randomly using Openssl RAND_bytes */
    unsigned int* AR;
    unsigned int* BR;
   
    
    AR = (unsigned int*)malloc(WORDLENGTH*sizeof(unsigned int));
    BR = (unsigned int*)malloc(WORDLENGTH*sizeof(unsigned int));
    for (size_t i=0; i < K/ (WORDLENGTH); i++) {
        RAND_pseudo_bytes(AR, 4*WORDLENGTH);
        RAND_pseudo_bytes(BR, 4*WORDLENGTH);
        if(compareArray(AR,BR,WORDLENGTH,WORDLENGTH)==-1){
            unsigned int* tempArr;
            tempArr=AR;
            AR=BR;
            BR=tempArr;
        }
        for(int j=0;j<WORDLENGTH;j++){
            A[i*WORDLENGTH+j]=AR[j];
            B[i*WORDLENGTH+j]=BR[j];
        }
    }
    
// for montgomery , ADDMOD  and SUBTRACTMOD test ensure A<PRIMEtest
    if(OPERATOR==MONTGOMERYMULTIPLICATION||OPERATOR==ADDMOD||OPERATOR==SUBTRACTMOD){
     for (size_t i=0; i < K; i+=WORDLENGTH){
    
       if (A[i]>=PRIME[0])    A[i]=PRIME[0]-1;
       if(B[i]>=PRIME[0])     B[i]=PRIME[0]-1;
     }
    }
    
    free(AR);
    free(BR);
    clock_gettime(CLOCK_MONOTONIC, &tend_init);
    
    
    /* Get platform/device information */
    ret = clGetPlatformIDs(4, platform_id, &ret_num_platforms);
  //  ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &ret_num_devices);
    ret = clGetDeviceIDs(platform_id[1], CL_DEVICE_TYPE_CPU, 1, &device_id, &ret_num_devices);
    
  //  printDeviceInfo(device_id);
    
    /* Create OpenCL Context */
    context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);
    
    /* Create command queue */
    command_queue = clCreateCommandQueue(context, device_id, 0, &ret);

    clock_gettime(CLOCK_MONOTONIC, &tend_createContext);
    
    /* Create buffer object */
    Amobj = clCreateBuffer(context, CL_MEM_READ_ONLY,  K*sizeof(unsigned int), NULL, &ret);
    Bmobj = clCreateBuffer(context, CL_MEM_READ_ONLY,  K*sizeof(unsigned int), NULL, &ret);
     if(OPERATOR==MULTIPLYOPRANDSCANNING||OPERATOR==MULTIPLYPRODUCTSCANNING) Cmobj = clCreateBuffer(context, CL_MEM_READ_WRITE, 2*K*sizeof(unsigned int), NULL, &ret);
    else Cmobj = clCreateBuffer(context, CL_MEM_READ_WRITE, K*sizeof(unsigned int), NULL, &ret);

    Omobj = clCreateBuffer(context, CL_MEM_READ_WRITE, 4*sizeof(int), NULL, &ret);
    Pmobj = clCreateBuffer(context, CL_MEM_READ_WRITE, WORDLENGTH*sizeof(unsigned int), NULL, &ret);

    /* Copy input data to memory buffer */
    ret = clEnqueueWriteBuffer(command_queue, Amobj, CL_TRUE, 0, K*sizeof(unsigned int), A, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, Bmobj, CL_TRUE, 0, K*sizeof(unsigned int), B, 0, NULL, NULL);
    OPERATOR_WORDSIZE_BITSLENGHT_MPRIME[0]=OPERATOR;
    OPERATOR_WORDSIZE_BITSLENGHT_MPRIME[1]=WORDSIZE;
    OPERATOR_WORDSIZE_BITSLENGHT_MPRIME[2]=BITSLENGTH;
    OPERATOR_WORDSIZE_BITSLENGHT_MPRIME[3]=MPRIME;

    ret = clEnqueueWriteBuffer(command_queue, Omobj, CL_TRUE, 0, 4*sizeof(int), OPERATOR_WORDSIZE_BITSLENGHT_MPRIME , 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, Pmobj, CL_TRUE, 0, WORDLENGTH*sizeof(unsigned int), PRIME, 0, NULL, NULL);

    clock_gettime(CLOCK_MONOTONIC, &tend_loadTomemory);
    
    /* Create kernel from source */
    program = clCreateProgramWithSource(context, 1, (const  char **)&source_str, (const size_t *)&source_size, &ret);
    ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);

    if (ret != CL_SUCCESS) {
        char buffer[10240];
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, NULL);
        fprintf(stderr, "CL Compilation failed:\n%s", buffer);
        abort();
    }
    
     clock_gettime(CLOCK_MONOTONIC, &tend_BuildProgram);
    /* Create task parallel OpenCL kernel */
    kernel = clCreateKernel(program, "mpaKernel", &ret);
    if (ret != CL_SUCCESS)
    {
        printf("Error: Failed to create kernel ! %s\n", getErrorString(ret));
        exit(1);
    }
    
    /* Set OpenCL kernel arguments */
    
        ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&Amobj);
        ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&Bmobj);
        ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&Cmobj);
        ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&Omobj);
        ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&Pmobj);
    
    if (ret != CL_SUCCESS)
    {
        printf("Error: Failed to set kernel arguments! %s\n", getErrorString(ret));
        exit(1);
    }
    clock_gettime(CLOCK_MONOTONIC, &tend_createKernel);
    /* Execute OpenCL kernel as task parallel */
    
        ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, global, local, 0, NULL, NULL);
        if (ret)
        {
            printf("Error: Failed to execute kernel %s!\n",getErrorString(ret));
            return EXIT_FAILURE;
        }
    
    clFinish(command_queue);

    clock_gettime(CLOCK_MONOTONIC, &tend_exec);
    /* Copy result to host */
   if(OPERATOR==MULTIPLYOPRANDSCANNING||OPERATOR==MULTIPLYPRODUCTSCANNING)  ret = clEnqueueReadBuffer(command_queue, Cmobj, CL_TRUE, 0, 2*K*sizeof(unsigned int), C, 0, NULL, NULL);
   else ret = clEnqueueReadBuffer(command_queue, Cmobj, CL_TRUE, 0, K*sizeof(unsigned int), C, 0, NULL, NULL);
    printf("clEnqueueReadBuffer for Cmobj  %s \n",getErrorString(ret));
    clFinish(command_queue);

    clock_gettime(CLOCK_MONOTONIC, &tend_redResults);
    //Display result 
   /* for (i=0; i < 256; i++) {
       
            printf("%u \n", C[K-i]);
       
      
    }*/
     printf("Entring Test for %s OPERATOR  using K=%zu and WORDLENGTH=%d and BITSLENGTH=%d \n",decode(OPERATOR),K,WORDLENGTH ,BITSLENGTH);
    if(testGpuResults(A,B,C,K,OPERATOR,DEBUG_MODE,PRIME,WORDLENGTH,bigPrime)==1)printf("%s executed %zu times Successefully \n",decode(OPERATOR),K, WORDLENGTH, bigPrime);
    clock_gettime(CLOCK_MONOTONIC, &tend_test);
double Initialization=0,CREATECONTEXT=0, LOADToMEMORY=0,BuildProgram=0,CREATEKERNEL=0,EXECUTION=0,READRESULTS=0,  CPUTIME=0,OPENCLOVRALLTime=0, SPEEDUP;

  Initialization =  ((double)tend_init.tv_sec + 1.0e-9*tend_init.tv_nsec) - 
           ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
 CREATECONTEXT =   ((double)tend_createContext.tv_sec + 1.0e-9*tend_createContext.tv_nsec) - 
           ((double)tend_init.tv_sec + 1.0e-9*tend_init.tv_nsec);
LOADToMEMORY =((double)tend_loadTomemory.tv_sec + 1.0e-9*tend_loadTomemory.tv_nsec) - 
           ((double)tend_createContext.tv_sec + 1.0e-9*tend_createContext.tv_nsec);
BuildProgram = ((double)tend_BuildProgram.tv_sec + 1.0e-9*tend_BuildProgram.tv_nsec) - 
           ((double)tend_loadTomemory.tv_sec + 1.0e-9*tend_loadTomemory.tv_nsec);
CREATEKERNEL = ((double)tend_createKernel.tv_sec + 1.0e-9*tend_createKernel.tv_nsec) - 
           ((double)tend_BuildProgram.tv_sec + 1.0e-9*tend_BuildProgram.tv_nsec),
EXECUTION = ((double)tend_exec.tv_sec + 1.0e-9*tend_exec.tv_nsec) - 
           ((double)tend_createKernel.tv_sec + 1.0e-9*tend_createKernel.tv_nsec),
READRESULTS = ((double)tend_redResults.tv_sec + 1.0e-9*tend_redResults.tv_nsec) - 
           ((double)tend_exec.tv_sec + 1.0e-9*tend_exec.tv_nsec),
CPUTIME =((double)tend_test.tv_sec + 1.0e-9*tend_test.tv_nsec) - 
           ((double)tend_redResults.tv_sec + 1.0e-9*tend_redResults.tv_nsec);
OPENCLOVRALLTime= LOADToMEMORY+ EXECUTION + READRESULTS;
      
    // dispaly timing in table 
    printf("----------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
 if(OPENCLOVRALLTime<CPUTIME){
    SPEEDUP = (CPUTIME/OPENCLOVRALLTime)*100;
    printf(  "Initialization  | CREATE CONTEXT  |  LOAD To MEMORY |Build Program src| CREATE KERNEL  |  EXECUTION      |      READ RESULTS    " RED_TERMINAL  "|    CPUTIME   " END_COLOR GREEN_TERMINAL "| OPENCLOVRALLTime |   SPEEDUP   |\n" END_COLOR);

}
else { printf(  "Initialization  | CREATE CONTEXT  |  LOAD To MEMORY |Build Program src| CREATE KERNEL  |  EXECUTION      |      READ RESULTS    " GREEN_TERMINAL  "|    CPUTIME   " END_COLOR RED_TERMINAL "| OPENCLOVRALLTime | CPU SPEEDUP |\n" END_COLOR);
    SPEEDUP = (OPENCLOVRALLTime/CPUTIME)*100;
    }
    printf(                  "    %.6f    |     %.6f    |     %.6f    |     %.6f    |    %.6f    |    %.6f     |      %.6f        |   %.6f    |    %.6f     |    %.2f %%    |\n", 
                         Initialization   ,  CREATECONTEXT    ,  LOADToMEMORY ,   BuildProgram  , CREATEKERNEL ,  EXECUTION,            READRESULTS      ,            CPUTIME,        OPENCLOVRALLTime, SPEEDUP  );
    printf("----------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    free(source_str);
    free(A);
    free(B);
    free(C);
    /* Finalization */
    ret = clFlush(command_queue);
    ret = clFinish(command_queue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseProgram(program);
    ret = clReleaseMemObject(Amobj);
    ret = clReleaseMemObject(Bmobj);
    ret = clReleaseMemObject(Cmobj);
    ret = clReleaseMemObject(Omobj);
    ret = clReleaseCommandQueue(command_queue);

     printf("clReleaseContext %s \n",getErrorString(ret));
    ret = clReleaseContext(context);

    
    return 0;
}

int testGpuResults(unsigned int* input1, unsigned int* input2, unsigned int* outputBytes, size_t K,  int OPERATOR, int DEBUG_MODE,  unsigned int* PRIME, int WORDLENGTH , mpz_t bigPrime){
int result=1;   
   // declare bigA and bigB outside loop 
   mpz_t bigA;
   mpz_t bigB;

   //Define TWOPOW_WL 
   mpz_t TWOPOW_WL;
   mpz_t RmoinsUn;
   mpz_init(TWOPOW_WL);
   mpz_init(RmoinsUn);
   const unsigned int* TWOPOW_WLStr="10000000000000000000000000000000000000000000000000000000000000000";
   mpz_set_str(TWOPOW_WL,TWOPOW_WLStr, 16);
   
   //mpz_out_str(stdout,16,TWOPOW_WL);
   switch(WORDLENGTH){
     case 512 : mpz_mul(TWOPOW_WL,TWOPOW_WL,TWOPOW_WL);
     case 256 : mpz_mul(TWOPOW_WL,TWOPOW_WL,TWOPOW_WL);
     case 128 : mpz_mul(TWOPOW_WL,TWOPOW_WL,TWOPOW_WL);
     case 64 : mpz_mul(TWOPOW_WL,TWOPOW_WL,TWOPOW_WL);
     case 32 : mpz_mul(TWOPOW_WL,TWOPOW_WL,TWOPOW_WL);
     case 16 : mpz_mul(TWOPOW_WL,TWOPOW_WL,TWOPOW_WL);
     
  

   }
   unsigned int* GMPBytes;
// for tsets to be removed
  /* mpz_sub_ui(TWOPOW_WL,TWOPOW_WL,1);
  mpz_out_str(stdout,16,TWOPOW_WL);
  printf("\n");
  */
 // mpz_out_str(stdout,16,bigPrime);
 // printf("\n");

    
   size_t* count;
   count = (size_t*) malloc(sizeof(size_t));

   unsigned int aBytes[WORDLENGTH];
   unsigned int bBytes[WORDLENGTH];
   unsigned int resultBytes[WORDLENGTH*2];
   
  if(OPERATOR==MULTIPLYOPRANDSCANNING||OPERATOR==MULTIPLYPRODUCTSCANNING) 
   
    GMPBytes=(unsigned int*)malloc(WORDLENGTH*2*sizeof(unsigned int));
    
  else 
    
    GMPBytes=(unsigned int*)malloc(WORDLENGTH*sizeof(unsigned int));
  
   


for(size_t i=0;i<K/WORDLENGTH;i++){
     if(OPERATOR==MULTIPLYOPRANDSCANNING||OPERATOR==MULTIPLYPRODUCTSCANNING) memcpy(resultBytes,&outputBytes[i*WORDLENGTH*2],2*WORDLENGTH*sizeof(unsigned int) );
   else memcpy(resultBytes,&outputBytes[i*WORDLENGTH],WORDLENGTH*sizeof(unsigned int) );
 
// copy the i th element

 //  unsigned long  a[20];
//mpz_import (z, 20, 1, sizeof(a[0]), 0, 0, a);
memcpy(aBytes,&input1[i*WORDLENGTH],sizeof(aBytes)  );
mpz_init(bigA);
mpz_import(bigA, WORDLENGTH, 1, sizeof(aBytes[0]), 0, 0, aBytes);

// copy the i th element
memcpy(bBytes,&input2[i*WORDLENGTH],WORDLENGTH*sizeof(unsigned int) );
mpz_init(bigB);
mpz_import(bigB, WORDLENGTH, 1, sizeof(bBytes[0]), 0, 0, bBytes);


/*
printArray(aBytes,WORDLENGTH,0);
mpz_out_str(stdout,16,bigA);printf("=BigA\n");
printArray(bBytes,WORDLENGTH,0);
mpz_out_str(stdout,16,bigB);printf("=BigB\n");
*/

if(OPERATOR==ADD){
// bigA=bigB+bigA
mpz_add(bigA,bigB,bigA); 
if(mpz_cmp(bigA,TWOPOW_WL)==1) mpz_sub(bigA,bigA,TWOPOW_WL);
// mpz_out_str(stdout,16,bigA);printf("=BigA*BigB\n");
}
if(OPERATOR==SUBTRACT){
// bigA=bigA-bigB%PRIME
 mpz_sub(bigA,bigA,bigB);
}
if(OPERATOR==MONTGOMERYMULTIPLICATION){
// R.modInverse(bigPrime).multiply(bigA).multiply(bigB).mod(bigPrime);
    //  int mpz_invert (mpz_t rop, const mpz_t op1, const mpz_t op2)
mpz_invert(RmoinsUn,TWOPOW_WL,bigPrime);
mpz_mul(bigA,RmoinsUn,bigA);
mpz_mul(bigA,bigA,bigB);
mpz_mod(bigA,bigA,bigPrime);
}
if(OPERATOR==SUBTRACTMOD){
// bigA=bigA-bigB%PRIME
 mpz_sub(bigA,bigA,bigB);
 mpz_mod(bigA,bigA,bigPrime);
}
if(OPERATOR==ADDMOD){
// bigA=(bigA+bigB)%PRIME
 mpz_add(bigA,bigA,bigB);
 mpz_mod(bigA,bigA,bigPrime);
}
if(OPERATOR==MULTIPLYOPRANDSCANNING||OPERATOR==MULTIPLYPRODUCTSCANNING){
// bigA=bigB+bigA
mpz_mul(bigA,bigB,bigA);

}
if((OPERATOR==MULTIPLYOPRANDSCANNING||OPERATOR==MULTIPLYPRODUCTSCANNING))
  mpz_export(GMPBytes, count, -1, 2*WORDLENGTH*sizeof(unsigned int), 0, 0, bigA);
else 
mpz_export(GMPBytes, count, -1, WORDLENGTH*sizeof(unsigned int), 0, 0, bigA);

//mpz_export((void*)GMPBytes, count, 1, sizeof(unsigned int), 1, 0, bigB);
// the bellow block should be uncommented in debug mod 
if((OPERATOR==MULTIPLYOPRANDSCANNING||OPERATOR==MULTIPLYPRODUCTSCANNING))
 { if(DEBUG_MODE!=0&&compareArray(resultBytes,GMPBytes,count[0],WORDLENGTH*2) != 0) {
    printf("64 mult Error at index %zu\n", i);
    result=0;
 printf("aBytes      = ");  printArray(aBytes,WORDLENGTH,0);
 printf("bBytes      = ");     printArray(bBytes,WORDLENGTH,0);
 printf("resultBytes = ");      printArray(resultBytes,WORDLENGTH*2,0);
 printf("GMPBytes    = ");     printArray(GMPBytes,WORDLENGTH*2,0);
 mpz_out_str(stdout,16,bigA);
  printf ("\n");
     break;
        } 
}
else if(DEBUG_MODE!=0&&compareArray(resultBytes,GMPBytes,count[0],WORDLENGTH) != 0) {
    printf("Error at index %zu\n", i);
    result=0;
 printf("PRIME      = ");      printArray(PRIME,WORDLENGTH,0);
 printf("aBytes      = ");     printArray(aBytes,WORDLENGTH,0);
 printf("bBytes      = ");     printArray(bBytes,WORDLENGTH,0);
 printf("resultBytes = ");     printArray(resultBytes,WORDLENGTH,0);
 printf("GMPBytes    = ");     printArray(GMPBytes,WORDLENGTH,0);
 mpz_out_str(stdout,16,bigA);
  printf ("\n");
     break;
        }

       
   
 
 
}
  mpz_clear(bigA);
  mpz_clear(bigB);
  mpz_clear(TWOPOW_WL);
  free(count);
  free(GMPBytes);
  
 
 return result;

}




