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
#define MULTIPLYOPERANDSCANNING 5
#define MULTIPLYPRODUCTSCANNING 6
#define MONTGOMERYMULTIPLICATION 7
#define ARITHMETICS 8

#define MAX_SOURCE_SIZE (0x100000)
typedef unsigned char WORDT;
#define  END_COLOR   "\x1b[0m"
#define  BLUE_TERMINAL    "\x1b[34m"
#define  RED_TERMINAL     "\x1b[31m"
#define  GREEN_TERMINAL   "\x1b[32m"

int compareArray(unsigned char*  a,unsigned char*  b,const int SIZE,const int WORDLINGTH) {
    int diff=WORDLINGTH-SIZE;
for (int i = 0; i < SIZE; i++){
if(a[i+diff]>b[i]) {
    return 1;}
if(a[i+diff]<b[i]) {
    return -1;
}
}
return 0;
}
char ansi[]= "\x1b[32m";
void printArray(unsigned char*  bytes,const int SIZE,const size_t ID) {

        char charAti[5];
        sprintf(charAti,"[%3ld", bytes[ID*SIZE]);
        printf("%s",charAti );
           for (int i = 1; i < SIZE; i++){

            sprintf(charAti,",%3ld", bytes[i+ID*SIZE]);
            printf("%s",charAti );
     }

        printf("]\n" );

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
        case MULTIPLYOPERANDSCANNING: return "MULTIPLYOPERANDSCANNING";
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

static void mpaToWords(const mpz_t z, WORDT *buf, int nwords, int wbits)
{
    mpz_t t;
    unsigned long mask = (wbits >= 32) ? 0xFFFFFFFFUL : ((1UL << wbits) - 1UL);
    int i;
    mpz_init_set(t, z);
    for (i = nwords - 1; i >= 0; i--) {
        buf[i] = (WORDT)(mpz_get_ui(t) & mask);
        mpz_tdiv_q_2exp(t, t, (mp_bitcnt_t)wbits);
    }
    mpz_clear(t);
}

static void mpaFromWords(mpz_t z, const WORDT *buf, int nwords, int wbits)
{
    int i;
    mpz_set_ui(z, 0);
    for (i = 0; i < nwords; i++) {
        mpz_mul_2exp(z, z, (mp_bitcnt_t)wbits);
        mpz_add_ui(z, z, (unsigned long)buf[i]);
    }
}

static unsigned long mpaMPrime(const mpz_t p, int wbits)
{
    mpz_t base, inv, mp;
    unsigned long r;
    mpz_inits(base, inv, mp, NULL);
    mpz_ui_pow_ui(base, 2, (unsigned long)wbits);
    if (mpz_invert(inv, p, base) == 0) {
        fprintf(stderr, "modulus is even; no Montgomery inverse exists\n");
        exit(EXIT_FAILURE);
    }
    mpz_sub(mp, base, inv);
    r = mpz_get_ui(mp);
    mpz_clears(base, inv, mp, NULL);
    return r;
}

static void mpaPickDevice(cl_device_id *outDev)
{
    cl_platform_id plats[16];
    cl_uint nplat = 0, nd = 0;
    cl_device_type order[3];
    const char *want = getenv("MPA_DEVICE_TYPE");
    int t, t0 = 0;
    cl_uint i;
    cl_device_id d;

    order[0] = CL_DEVICE_TYPE_GPU;
    order[1] = CL_DEVICE_TYPE_ACCELERATOR;
    order[2] = CL_DEVICE_TYPE_CPU;

    if (clGetPlatformIDs(16, plats, &nplat) != CL_SUCCESS || nplat == 0) {
        fprintf(stderr, "no OpenCL platform found\n");
        exit(EXIT_FAILURE);
    }
    if (want && !strcmp(want, "cpu")) t0 = 2;

    for (t = t0; t < 3; t++)
        for (i = 0; i < nplat; i++)
            if (clGetDeviceIDs(plats[i], order[t], 1, &d, &nd) == CL_SUCCESS && nd > 0) {
                *outDev = d;
                return;
            }
    fprintf(stderr, "no OpenCL device found\n");
    exit(EXIT_FAILURE);
}

static int testGpuResults(const WORDT *input1, const WORDT *input2,
                          const WORDT *outputBytes, size_t K, int OPERATOR,
                          int DEBUG_MODE, const WORDT *PRIME, int WORDLENGTH,
                          const mpz_t bigPrime, int wbits);

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

           errno = 0;
           iterations = strtoull(str, &endptr, base);

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
       unsigned long MPRIME;

    WORDT* PRIME;
    PRIME= (WORDT*)malloc(WORDLENGTH*sizeof(WORDT));
    if (!PRIME) { fprintf(stderr, "out of memory\n"); exit(EXIT_FAILURE); }
     mpz_t bigPrime;
     mpz_init(bigPrime);
    {
    const char *primeStr = NULL;
    switch(BITSLENGTH){
        case 256 :
            primeStr = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";
            break;
        case 512 :
            primeStr = "AADD9DB8DBE9C48B3FD4E6AE33C9FC07CB308DB3B3C9D20ED6639CCA70330871"
                       "7D4D9B009BC66842AECDA12AE6A380E62881FF2F2D82C68528AA6056583A48F3";
            break;
        case 1024 :
            primeStr = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                       "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                       "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                       "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF97";
            break;
        case 2048 :
            primeStr = "E53DF5FC3F650D066875837012A4E7BEA863C65CB592D9C36942CF69CBC6DD4F"
                       "D804E19CCF2696C9BEBCF18742FA5FB091CBDE1782E8291009464913ECE19745"
                       "7800EA6E43B0E2A64615D182B6DE150479C58D1C7C702D47EA3031B379CA13A2"
                       "048C964E1D1E8D4CD3815D0895BF31E53271D4607E16461B77FB26100915D679"
                       "9060203EDEBFEA9495A5A8E7CED68FC9DB2D47CE7992461BA78174608AD0BBE3"
                       "F5E63EC6C960564430CBD2E6E587D08EE12F94B5B99DFFB12C6727A25E800DAC"
                       "6CD8DE77A5BBC93B36E444B070888CB5ADD991870466968A6E9A23C2EE0A1D67"
                       "1C9B601081A44AA6A58D4DC76686EF15FCE1C9AEB4033395A9B24BE1AA1929BB";
            break;
        default :
            fprintf(stderr, "unsupported BITSLENGTH %d (expected 256, 512, 1024 or 2048)\n",
                    BITSLENGTH);
            exit(EXIT_FAILURE);
    }
    if (mpz_set_str(bigPrime, primeStr, 16) != 0) {
        fprintf(stderr, "malformed modulus literal\n"); exit(EXIT_FAILURE);
    }
    if ((int)mpz_sizeinbase(bigPrime, 2) != BITSLENGTH) {
        fprintf(stderr, "modulus literal is %d bits, expected %d\n",
                (int)mpz_sizeinbase(bigPrime, 2), BITSLENGTH);
        exit(EXIT_FAILURE);
    }
    mpaToWords(bigPrime, PRIME, WORDLENGTH, WORDSIZE);
    MPRIME = mpaMPrime(bigPrime, WORDSIZE);
    }

    struct timespec tstart={0,0}, tend_init={0,0} , tend_createContext={0,0},
    tend_loadTomemory={0,0},tend_BuildProgram={0,0}, tend_createKernel={0,0}, tend_exec={0,0}, tend_redResults={0,0}, tend_test={0,0};
    clock_gettime(CLOCK_MONOTONIC, &tstart);

    cl_platform_id platform_id = NULL;
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
   const size_t global[]={K/WORDLENGTH};
   const size_t *local = NULL;
    int i, j;
    unsigned char* A;
    unsigned char* B;
    unsigned char* C;
    int* OPERATOR_WORDSIZE_BITSLENGHT_MPRIME;

    A = (unsigned char*)malloc(K*sizeof(unsigned char));
    B = (unsigned char*)malloc(K*sizeof(unsigned char));
    if(OPERATOR==MULTIPLYOPERANDSCANNING||OPERATOR==MULTIPLYPRODUCTSCANNING) C = (unsigned char*)malloc(2*K*sizeof(unsigned char));
    else C = (unsigned char*)malloc(K*sizeof(unsigned char));
    OPERATOR_WORDSIZE_BITSLENGHT_MPRIME = (int*)malloc(4*sizeof(int));

    FILE *fp;
    const char fileName[] = "mpaKernels_8bits.cl";
    size_t source_size;
     char *source_str;

    fp = fopen(fileName, "rb");
    if (!fp) {
        fprintf(stderr, "Failed to load kernel.\n");
        exit(1);
    }
    source_str = ( char *)malloc(MAX_SOURCE_SIZE);
    source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
    fclose(fp);

    unsigned char* AR;
    unsigned char* BR;

    AR = (unsigned char*)malloc(WORDLENGTH*sizeof(unsigned char));
    BR = (unsigned char*)malloc(WORDLENGTH*sizeof(unsigned char));
    for (size_t i=0; i < K/ (WORDLENGTH); i++) {
        if (RAND_bytes((unsigned char *)AR, (int)(WORDLENGTH)) != 1) { fprintf(stderr, "RAND_bytes failed\n"); exit(EXIT_FAILURE); }
        if (RAND_bytes((unsigned char *)BR, (int)(WORDLENGTH)) != 1) { fprintf(stderr, "RAND_bytes failed\n"); exit(EXIT_FAILURE); }
        if(compareArray(AR,BR,WORDLENGTH,WORDLENGTH)==-1){
            unsigned char* tempArr;
            tempArr=AR;
            AR=BR;
            BR=tempArr;
        }
        for(int j=0;j<WORDLENGTH;j++){
            A[i*WORDLENGTH+j]=AR[j];
            B[i*WORDLENGTH+j]=BR[j];
        }
    }

    if(OPERATOR==MONTGOMERYMULTIPLICATION||OPERATOR==ADDMOD||OPERATOR==SUBTRACTMOD){
     for (size_t i=0; i < K; i+=WORDLENGTH){

       if (A[i]>=PRIME[0])    A[i]=PRIME[0]-1;
       if(B[i]>=PRIME[0])     B[i]=PRIME[0]-1;
     }
    }

    free(AR);
    free(BR);
    clock_gettime(CLOCK_MONOTONIC, &tend_init);

    ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    mpaPickDevice(&device_id);

    context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);

    command_queue = clCreateCommandQueue(context, device_id, 0, &ret);

    clock_gettime(CLOCK_MONOTONIC, &tend_createContext);

    Amobj = clCreateBuffer(context, CL_MEM_READ_ONLY,  K*sizeof(unsigned char), NULL, &ret);
    Bmobj = clCreateBuffer(context, CL_MEM_READ_ONLY,  K*sizeof(unsigned char), NULL, &ret);
     if(OPERATOR==MULTIPLYOPERANDSCANNING||OPERATOR==MULTIPLYPRODUCTSCANNING) Cmobj = clCreateBuffer(context, CL_MEM_READ_WRITE, 2*K*sizeof(unsigned char), NULL, &ret);
    else Cmobj = clCreateBuffer(context, CL_MEM_READ_WRITE, K*sizeof(unsigned char), NULL, &ret);

    Omobj = clCreateBuffer(context, CL_MEM_READ_WRITE, 4*sizeof(int), NULL, &ret);
    Pmobj = clCreateBuffer(context, CL_MEM_READ_WRITE, WORDLENGTH*sizeof(unsigned char), NULL, &ret);

    ret = clEnqueueWriteBuffer(command_queue, Amobj, CL_TRUE, 0, K*sizeof(unsigned char), A, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, Bmobj, CL_TRUE, 0, K*sizeof(unsigned char), B, 0, NULL, NULL);
    OPERATOR_WORDSIZE_BITSLENGHT_MPRIME[0]=OPERATOR;
    OPERATOR_WORDSIZE_BITSLENGHT_MPRIME[1]=WORDSIZE;
    OPERATOR_WORDSIZE_BITSLENGHT_MPRIME[2]=BITSLENGTH;
    OPERATOR_WORDSIZE_BITSLENGHT_MPRIME[3]=(int)(unsigned int)MPRIME;

    ret = clEnqueueWriteBuffer(command_queue, Omobj, CL_TRUE, 0, 4*sizeof(int), OPERATOR_WORDSIZE_BITSLENGHT_MPRIME , 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, Pmobj, CL_TRUE, 0, WORDLENGTH*sizeof(unsigned char), PRIME, 0, NULL, NULL);

    clock_gettime(CLOCK_MONOTONIC, &tend_loadTomemory);

    program = clCreateProgramWithSource(context, 1, (const  char **)&source_str, (const size_t *)&source_size, &ret);
    char buildOpts[128];
    snprintf(buildOpts, sizeof(buildOpts), "-I%s -DWORDLENGTH_T=%d",
             getenv("MPA_KERNEL_DIR") ? getenv("MPA_KERNEL_DIR") : ".", WORDLENGTH);
    ret = clBuildProgram(program, 1, &device_id, buildOpts, NULL, NULL);

    if (ret != CL_SUCCESS) {
        char buffer[10240];
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, NULL);
        fprintf(stderr, "CL Compilation failed:\n%s", buffer);
        abort();
    }

     clock_gettime(CLOCK_MONOTONIC, &tend_BuildProgram);
    kernel = clCreateKernel(program, "mpaKernel", &ret);
    if (ret != CL_SUCCESS)
    {
        printf("Error: Failed to create kernel ! %s\n", getErrorString(ret));
        exit(1);
    }

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

        ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, global, local, 0, NULL, NULL);
        if (ret)
        {
            printf("Error: Failed to execute kernel %s!\n",getErrorString(ret));
            return EXIT_FAILURE;
        }

    clFinish(command_queue);

    clock_gettime(CLOCK_MONOTONIC, &tend_exec);
   if(OPERATOR==MULTIPLYOPERANDSCANNING||OPERATOR==MULTIPLYPRODUCTSCANNING)  ret = clEnqueueReadBuffer(command_queue, Cmobj, CL_TRUE, 0, 2*K*sizeof(unsigned char), C, 0, NULL, NULL);
   else ret = clEnqueueReadBuffer(command_queue, Cmobj, CL_TRUE, 0, K*sizeof(unsigned char), C, 0, NULL, NULL);
    printf("clEnqueueReadBuffer for Cmobj  %s \n",getErrorString(ret));
    clFinish(command_queue);

    clock_gettime(CLOCK_MONOTONIC, &tend_redResults);

     printf("Entring Test for %s OPERATOR  using K=%zu and WORDLENGTH=%d and BITSLENGTH=%d \n",decode(OPERATOR),K,WORDLENGTH ,BITSLENGTH);
    int verified = testGpuResults(A,B,C,K,OPERATOR,DEBUG_MODE,PRIME,WORDLENGTH,bigPrime,WORDSIZE);
    if (verified)
        printf(GREEN_TERMINAL "%s verified against GMP for all %zu items" END_COLOR "\n",
               decode(OPERATOR), K/(size_t)WORDLENGTH);
    else
        printf(RED_TERMINAL "%s FAILED verification" END_COLOR "\n", decode(OPERATOR));
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
    ret = clFlush(command_queue);
    ret = clFinish(command_queue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseProgram(program);
    ret = clReleaseMemObject(Amobj);
    ret = clReleaseMemObject(Bmobj);
    ret = clReleaseMemObject(Cmobj);
    ret = clReleaseMemObject(Omobj);
    ret = clReleaseMemObject(Pmobj);
    ret = clReleaseCommandQueue(command_queue);

    ret = clReleaseContext(context);

    free(PRIME);

    free(OPERATOR_WORDSIZE_BITSLENGHT_MPRIME);
    mpz_clear(bigPrime);
    return verified ? 0 : 1;
}

static int testGpuResults(const WORDT *input1, const WORDT *input2,
                          const WORDT *outputBytes, size_t K, int OPERATOR,
                          int DEBUG_MODE, const WORDT *PRIME, int WORDLENGTH,
                          const mpz_t bigPrime, int wbits)
{
    const int outWords = (OPERATOR == MULTIPLYOPERANDSCANNING ||
                          OPERATOR == MULTIPLYPRODUCTSCANNING)
                         ? 2 * WORDLENGTH : WORDLENGTH;
    const size_t items = K / (size_t)WORDLENGTH;
    size_t i;
    long bad = 0;
    int rc = 1;

    mpz_t a, b, want, lim, R, Rinv;
    WORDT *expect = (WORDT *)malloc((size_t)outWords * sizeof(WORDT));
    if (!expect) { fprintf(stderr, "out of memory\n"); exit(EXIT_FAILURE); }

    mpz_inits(a, b, want, lim, R, Rinv, NULL);
    mpz_ui_pow_ui(lim, 2, (unsigned long)(wbits * WORDLENGTH));
    mpz_set(R, lim);
    if (OPERATOR == MONTGOMERYMULTIPLICATION) mpz_invert(Rinv, R, bigPrime);

    for (i = 0; i < items; i++) {
        int w, ok = 1;
        mpaFromWords(a, &input1[i * (size_t)WORDLENGTH], WORDLENGTH, wbits);
        mpaFromWords(b, &input2[i * (size_t)WORDLENGTH], WORDLENGTH, wbits);

        switch (OPERATOR) {
        case ADD:         mpz_add(want, a, b); mpz_mod(want, want, lim); break;
        case SUBTRACT:    mpz_sub(want, a, b); mpz_mod(want, want, lim); break;
        case ADDMOD:      mpz_add(want, a, b); mpz_mod(want, want, bigPrime); break;
        case SUBTRACTMOD: mpz_sub(want, a, b); mpz_mod(want, want, bigPrime); break;
        case MULTIPLYOPERANDSCANNING:
        case MULTIPLYPRODUCTSCANNING:
                          mpz_mul(want, a, b); break;
        case MONTGOMERYMULTIPLICATION:
                          mpz_mul(want, a, b);
                          mpz_mul(want, want, Rinv);
                          mpz_mod(want, want, bigPrime);
                          break;
        default:
            fprintf(stderr, "no reference for operator %d\n", OPERATOR);
            rc = 0; goto done;
        }

        mpaToWords(want, expect, outWords, wbits);
        for (w = 0; w < outWords; w++)
            if (outputBytes[i * (size_t)outWords + w] != expect[w]) { ok = 0; break; }

        if (!ok) {
            bad++;
            rc = 0;
            if (DEBUG_MODE && bad <= 3) {
                printf(RED_TERMINAL "mismatch at item %zu" END_COLOR "\n", i);
                gmp_printf("  a    = %Zx\n  b    = %Zx\n  want = %Zx\n", a, b, want);
            }
        }
    }

done:
    if (bad)
        printf(RED_TERMINAL "%s: %ld of %zu results WRONG" END_COLOR "\n",
               decode(OPERATOR), bad, items);
    mpz_clears(a, b, want, lim, R, Rinv, NULL);
    free(expect);
    return rc;
}
