/*----------------------------------------------------------------------
* Copyright 2022 Verasonics, Inc
*
* Permission is hereby granted, free of charge, to any
* person obtaining a copy of this software and associated
* documentation files (the "Software"), to deal in the
* Software without restriction, including without
* limitation the rights to use, copy, modify, merge,
* publish, distribute, sublicense, and/or sell copies of
* the Software, and to permit persons to whom the Software
* is furnished to do so, subject to the following
* conditions:
*
* The above copyright notice and this permission notice
* shall be included in all copies or substantial portions
* of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
* ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
* TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
* PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
* SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
* OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
* IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
**************************************************************************
* externalDoppler_GPU_MEXCUDA.cu
*   Example showing CDI processing using the GPU with CUDA code
*   Requires:
*       - NVidia GPU installed
*       - NVidia Driver for GPU installed
*       - CUDA 11.0 SDK to be installed in default location
*       - Parallel Computing Toolbox to be installed if mexcuda compilation
*         method is used
*
**************************************************************************/
#define UNIFIEDMEMORY 0  // Use unified memory (1) or not (0)
#define DEBUGGING 0      // Show debugging info in matlab terminal
#define SHOWTIMING 0     // Show processing timing info in matlab terminal

#define PI 3.14159
#define VELOCITYSCALE (-(253)/(2*PI)) //scale the phase measurement between -pi/pi and color 0 to 255

// includes
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// includes
#include <cuda_runtime.h>
#include <helper_cuda.h>  //helper functions found in CUDA sample code folder
#include <helper_functions.h>//helper functions found in CUDA sample code folder

// Structure definition for internal mxArray.  This definition allows setting the
//    number_array.pdata field manually to point to memory data that has no header.
typedef struct {
    void *reserved;
    int reserved1[2];
    void *reserved2;
    size_t number_of_dims;
    unsigned int reserved3;
    struct {
        unsigned int flag0 : 1;
        unsigned int flag1 : 1;
        unsigned int flag2 : 1;
        unsigned int flag3 : 1;
        unsigned int flag4 : 1;
        unsigned int flag5 : 1;
        unsigned int flag6 : 1;
        unsigned int flag7 : 1;
        unsigned int flag7a: 1;
        unsigned int flag8 : 1;
        unsigned int flag9 : 1;
        unsigned int flag10 : 1;
        unsigned int flag11 : 4;
        unsigned int flag12 : 8;
        unsigned int flag13 : 8;
    } flags;
    size_t reserved4[2];
    union {
        struct {
            void *pdata;
            void *pimag_data;
            void *reserved5;
            size_t reserved6[3];
        } number_array;
        struct {
            mxArray **pdata;
            char  *field_names;
            void  *dummy1;
            void  *dummy2;
            int   dummy3;
            int   nfields;
            } struct_array;
    } data;
} Internal_mxArray;

/* Variables */
static double *filterMatrixI;
static double *d_filterMatrixI=NULL;
static double *d_acI, *d_acQ;
static double *d_dataI, *d_dataQ;
static double *d_dataWFI, *d_dataWFQ;
static double *d_velocity;
static double *d_powerSq;
static double *velocity;
static double *powerSq;
static mxArray *mxCmdPrhs[3];

void Cleanup();

void Initialize(int numPixels, int Npri, int nrhs, mxArray const *prhs[], size_t memSize);

static __global__ void multiplyMatrices(double *mat1I, //matrix 1
                      double *mat1Q, // matrix 1
                      double *mat2I, // matrix 2
                      double *resI,  // result matrix
                      double *resQ,  // result matrix
                      int nrmat1,    // number of rows for matrix 1
                      int ncmat1,    // number of cols for matrix 1
                      int nrcmat2);  // number of rows & cols for matrix 2 (square)

static __global__ void lagoneautocorr(double *IQDataWFI,
                    double *IQDataWFQ,
                    double *autocorrI,
                    double *autocorrQ,
                    int numPixels,
                    int nPri);

void dispMatrix(double *matrixI, double *matrixQ,int nr, int nc);

static __global__ void calculatePowerAndVelocity(double *powerSq, double *velocity, double *acI, double *acQ, mwSize numPixels);


/* Main Function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    // Declare Variables
    const mwSize *dataSize;
    mwSize dataNumDims, numPixels, Npri;
    size_t memSize;
    double *IDataIn, *QDataIn;
    const mxArray *threshParams;
    mxArray *outImgData; // Pointer to output data (DOUBLE)

    //Check input parameters to function to prevent crashing
    if ( nrhs < 2)
    {
        mexErrMsgIdAndTxt("MyMex:gputest:nrhs", "Minimum of two inputs required, I & Q.");
    }
    if (nlhs!=1)
    {
        mexErrMsgIdAndTxt("MyMex:gputest:nlhs", "one output is required");
    }

    /* Get data sizes */
    dataNumDims = mxGetNumberOfDimensions(prhs[0]);
    dataSize = mxGetDimensions(prhs[0]);

    if (dataNumDims == 4)
    {
        numPixels = dataSize[0]*dataSize[1]*dataSize[2];
        Npri = dataSize[3];
    }
    else if (dataNumDims == 3)
    {
        numPixels = dataSize[0]*dataSize[1];
        Npri = dataSize[2];
    }

    // Convenient pointer to input data (DOUBLE)
    IDataIn = mxGetDoubles(prhs[0]);
    QDataIn = mxGetDoubles(prhs[1]);

    /* Calculate memory size */
    memSize = numPixels*Npri*sizeof(double);

    // if datasize has changed or wasnt known before, reinitialize filters
    if (d_filterMatrixI == NULL)
    {
        Initialize(numPixels, Npri, nrhs, prhs, memSize);
    }

    //////////// Processing ////////////////////
    /* Setup timing events */
    cudaEvent_t start, stop;
    float memcpyH2DTime_ms, memcpyD2HTime_ms, wallFilterTime_ms, autocorrTime_ms, powerThresTime_ms, calcPowerVel_ms;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    /* copy data to device */
    cudaEventRecord(start,0);
    checkCudaErrors(cudaMemcpyAsync(d_dataI, IDataIn, memSize, cudaMemcpyHostToDevice,0));
    checkCudaErrors(cudaMemcpyAsync(d_dataQ, QDataIn, memSize, cudaMemcpyHostToDevice,0));
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&memcpyH2DTime_ms, start, stop);

    /* wallfiltering : multiply data by filterMatrix*/
    cudaEventRecord(start,0);
    multiplyMatrices<<<(numPixels+255)/256, 256>>>(
            d_dataI, //matrix 1
            d_dataQ, //matrix 1
            d_filterMatrixI, //matrix 2
            d_dataWFI,  //result matrix
            d_dataWFQ,  //result matrix
            numPixels,   //number of rows for matrix 1
            Npri,   //number of cols for matrix 1
            Npri);  //number of rows & cols for matrix 2 (square)
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&wallFilterTime_ms, start, stop);

    #if DEBUGGING
        mexPrintf("----------- Data wall filtered --------------\n");
        dispMatrix(d_dataWFI, d_dataWFQ, numPixels, Npri);        //for debugging
    #endif


    /* lag-1 autocorrelation */
    cudaEventRecord(start,0);
    lagoneautocorr<<<(numPixels+255)/256, 256>>>(d_dataWFI, d_dataWFQ, d_acI,  d_acQ, numPixels, Npri);
    //The first argument in the execution configuration specifies the number of thread blocks in the grid,
    //and the second specifies the number of threads in a thread block.
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&autocorrTime_ms, start, stop);
    #if DEBUGGING
        mexPrintf("----------- Data autocorrelated --------------\n");
        dispMatrix(d_acI, d_acQ, numPixels, 1);        //for debugging
    #endif

     /* Threshold Data */
     cudaEventRecord(start,0);
     calculatePowerAndVelocity<<<(numPixels+255)/256, 256>>>(d_powerSq, d_velocity, d_acI, d_acQ, numPixels);
     cudaEventRecord(stop,0);
     cudaEventSynchronize(stop);
     cudaEventElapsedTime(&calcPowerVel_ms, start, stop);

     cudaEventRecord(start,0);
     checkCudaErrors(cudaMemcpyAsync(velocity, d_velocity, numPixels*sizeof(double), cudaMemcpyDeviceToHost,0));
     checkCudaErrors(cudaMemcpyAsync(powerSq, d_powerSq, numPixels*sizeof(double), cudaMemcpyDeviceToHost,0));
     cudaEventRecord(stop,0);
     cudaEventSynchronize(stop);
     cudaEventElapsedTime(&memcpyD2HTime_ms, start, stop);

     cudaEventRecord(start,0);
     threshParams = mexGetVariablePtr("base", "threshParams");
     outImgData = mxCreateNumericMatrix(numPixels, 1, mxDOUBLE_CLASS, mxREAL); // cannot preallocate because of Matlab memory
     mxCmdPrhs[2] = (mxArray *) threshParams;
     mexCallMATLAB(1, &outImgData, 3, mxCmdPrhs, "adaptiveThreshold");
     mxSetDimensions(outImgData, dataSize, dataNumDims-1); // reshape the dimensions for the output matrix
     cudaEventRecord(stop,0);
     cudaEventSynchronize(stop);
     cudaEventElapsedTime(&powerThresTime_ms, start, stop);
     plhs[0] = outImgData;


    #if SHOWTIMING
    mexPrintf("Memcpy Host-to-Device: %f ms\n", memcpyH2DTime_ms);
    mexPrintf("Wall filter: %f ms\n", wallFilterTime_ms);
    mexPrintf("Autocorrelation: %f ms\n", autocorrTime_ms);
    mexPrintf("Calc Power & Velocity: %f ms\n", calcPowerVel_ms);
    mexPrintf("Memcpy Device-to-Host: %f ms\n", memcpyD2HTime_ms);
    mexPrintf("Calculate Power & Threshold: %f ms\n", powerThresTime_ms);
    #endif
}

/* Cleanup Routine */
void cleanup()
{
    mexPrintf("Cleaning allocated memory for ExternalDoppler_GPU_MEXCUDA.mex!\n");
    checkCudaErrors(cudaFree(d_filterMatrixI));
    checkCudaErrors(cudaFree(d_acI));
    checkCudaErrors(cudaFree(d_acQ));
    checkCudaErrors(cudaFree(d_dataI));
    checkCudaErrors(cudaFree(d_dataQ));
    checkCudaErrors(cudaFree(d_dataWFI));
    checkCudaErrors(cudaFree(d_dataWFQ));
    checkCudaErrors(cudaFree(d_velocity));
    checkCudaErrors(cudaFree(d_powerSq));
}

void Initialize(int numPixels, int Npri, int nrhs, mxArray const *prhs[], size_t memSize)
{
    int i, ord;
            //int wfNumDims;
    double *wfIn;
    const mwSize *wfDims;
    char polystr[80];
    const mxArray *wallFilterWorkspace;

    filterMatrixI = (double*)mxCalloc(Npri*Npri, sizeof(double));
    ord = (int)ceil(Npri/8.0);
    mexPrintf("Initialize Wall Filter of order: %li\n", ord);

    // Allocate memory on device
    #if UNIFIEDMEMORY
    checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_filterMatrixI), Npri*Npri*sizeof(double)));
    checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_dataI), memSize));
    checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_dataQ), memSize));
    checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_dataWFI), memSize));
    checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_dataWFQ), memSize));
    checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_acI), numPixels*sizeof(double)));
    checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_acQ), numPixels*sizeof(double)));
    checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_velocity), numPixels*sizeof(double)));
    checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_powerSq), numPixels*sizeof(double)));
    #else
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_filterMatrixI), Npri*Npri*sizeof(double)));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_dataI), memSize));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_dataQ), memSize));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_dataWFI), memSize));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_dataWFQ), memSize));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_acI), numPixels*sizeof(double)));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_acQ), numPixels*sizeof(double)));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_velocity), numPixels*sizeof(double)));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_powerSq), numPixels*sizeof(double)));
    #endif
    checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&velocity), numPixels*sizeof(double)));
    checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&powerSq), numPixels*sizeof(double)));
    mexMakeMemoryPersistent(d_filterMatrixI);
    mexMakeMemoryPersistent(d_dataI);
    mexMakeMemoryPersistent(d_dataQ);
    mexMakeMemoryPersistent(d_dataWFI);
    mexMakeMemoryPersistent(d_dataWFQ);
    mexMakeMemoryPersistent(d_acI);
    mexMakeMemoryPersistent(d_acQ);
    mexMakeMemoryPersistent(d_velocity);
    mexMakeMemoryPersistent(d_powerSq);


    // allocate memory on host for thresholding
    mxCmdPrhs[0] = mxCreateNumericMatrix(numPixels, 1, mxDOUBLE_CLASS, mxREAL); //consider preallocating
    mxCmdPrhs[1] = mxCreateNumericMatrix(numPixels, 1, mxDOUBLE_CLASS, mxREAL); //consider preallocating
    mxFree(mxGetData(mxCmdPrhs[0])); // free the memory allocated by mxCreateNumericArray
    mxFree(mxGetData(mxCmdPrhs[1])); // free the memory allocated by mxCreateNumericArray
    ((Internal_mxArray*)(mxCmdPrhs[0]))->data.number_array.pdata = velocity;
    ((Internal_mxArray*)(mxCmdPrhs[1]))->data.number_array.pdata = powerSq;
    mexMakeArrayPersistent(mxCmdPrhs[0]);
    mexMakeArrayPersistent(mxCmdPrhs[1]);

    mexAtExit(cleanup);

    /* -- Wallfilter is an input -- */
    if (nrhs == 3)  //wallfilter has been provided, check size
    {
        wfIn = mxGetDoubles(prhs[2]); // Wallfilter
        wfDims = mxGetDimensions(prhs[2]);

        if (wfDims[0] != Npri)//check size of wallfilter
        {
            mexErrMsgIdAndTxt("MyMex:gputest:wallfilter", "both dimensions of the wallfilter(%i,%i) must equal the PRI(%i) of the ensemble", wfDims[0], wfDims[1],Npri);
        }
        for (i=0; i<(Npri*Npri); i++)
        {
            d_filterMatrixI[i] = (float)wfIn[i];
        }
    }
    else //generate walfilter in Matlab
    {
        // Initialize Wallfilter data using the Verasonics polyfilter.p function

        sprintf(polystr, "wallFilterMatrix = polyfilter(%i, %i);", Npri, ord);
        if (!mexEvalString(polystr))
        {
            mexPrintf("Using wall filter created in Matlab\n");
            wallFilterWorkspace = mexGetVariable("base","wallFilterMatrix");
            wfIn = mxGetDoubles(wallFilterWorkspace);
            for (i=0; i<(Npri*Npri); i++)
            {
                filterMatrixI[i] = (float)wfIn[i];
            }
        }
    }
    checkCudaErrors(cudaMemcpyAsync(d_filterMatrixI, filterMatrixI, Npri*Npri*sizeof(double), cudaMemcpyHostToDevice,0));
    cudaDeviceSynchronize();

    #if DEBUGGING
        mexPrintf("----------- WF --------------\n");
        dispMatrix(d_filterMatrixI,d_filterMatrixI, Npri, Npri);        //for debugging
    #endif
}

/*********************************************
 *  multiplyMatrices()
 *  function to multiply two Complex Matrices
 *
 **********************************************/

static __global__ void multiplyMatrices(
        double *mat1I, //matrix 1
        double *mat1Q, //matrix 1
        double *mat2I, //matrix 2
        double *resI,  //result matrix
        double *resQ,  //result matrix
        int nrmat1,   //number of rows for matrix 1
        int ncmat1,   //number of cols for matrix 1
        int nrcmat2)  //number of rows & cols for matrix 2 (square)
{

    #if DEBUGGING
        mexPrintf("nrmat1: %i ncmat1: %i nrcmat2: %i\n", nrmat1, ncmat1, nrcmat2);
    #endif

    int i, j, k;
    int a, b, c;
    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
    int n;
    for (n = threadID; n < (nrmat1*ncmat1); n += numThreads)
    {
        i = floor((float)n/(float)nrmat1);
        j = n % nrmat1;//mod(n,nrmat1);

        resI[n] = 0;
        resQ[n] = 0;

        c = i*nrcmat2;
        for (k = 0; k < nrcmat2; k++)
        {
            a = k*nrmat1 + j;
            b = c + k;
            // complex multiply = (x + yi)(u + vi) = (xu – yv) + //real
            //                                       (xv + yu)i  //imag
            resI[n] += mat1I[a]*mat2I[b];
                     //- mat1Q[a]*mat2Q[b]; //can leave off second part because filter is only real

            resQ[n] += mat1Q[a]*mat2I[b];
                     //+ mat1I[k*nrmat1 + j]*mat2Q[i*nrcmat2+k] + can leave off this part because filter is only real
        }
    }
}

/*****************************
 *  lagoneautocorr()
 *  Lag 1 autocorrelation
 *
 ******************************/
static __global__ void lagoneautocorr(
        double *IQDataWFI,
        double *IQDataWFQ,
        double *autocorrI,
        double *autocorrQ,
        int numPixels,
        int nPri)
{
    int a, b;
    double N = (double)(nPri-1);

    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

    for (int i = threadID; i < numPixels; i += numThreads)

    {
        autocorrI[i] = 0;
        autocorrQ[i] = 0;
        for (int j=0;j<(nPri-1);j++)
        {
            a = i+j*numPixels;
            b = a+numPixels;
            autocorrI[i] +=  (IQDataWFI[a])*( IQDataWFI[b]) +
                              (IQDataWFQ[a])*( IQDataWFQ[b]);

            autocorrQ[i] +=  (IQDataWFI[a])*( IQDataWFQ[b]) -
                              (IQDataWFI[b])*( IQDataWFQ[a]);
        }
        autocorrI[i]=autocorrI[i]/N;//calculate mean
        autocorrQ[i]=autocorrQ[i]/N;//calculate mean
    }
}

/**************************************************
 *  dispMatrix()
 *     display the data matrix for debugging
 *
 ***************************************************/
void dispMatrix(double *matrixI, double *matrixQ, int nr, int nc)
{
    int i, j;

    for (i=0;i<nr;i++)
    {
        for (j=0;j<nc;j++)
        {
            mexPrintf("%.1f + %.1f, \t", matrixI[nr*j+i], matrixQ[nr*j+i]);
        }
        mexPrintf("\n");
    }
}

/******************************
 *  calculatePowerAndVelocity()
 *
 ******************************/
static __global__ void calculatePowerAndVelocity(double *powerSq, double *velocity, double *acI, double *acQ, mwSize numPixels)
{
    int i;
    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

    for (i = threadID; i < numPixels; i += numThreads)

    {
        powerSq[i] = pow((acI[i]*acI[i] + acQ[i]*acQ[i]), 0.25);
        velocity[i] = VELOCITYSCALE * atan2(acQ[i], acI[i]);
    }
}