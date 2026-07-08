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
*
**************************************************************************
* externalFIR_GPU_MEXCUDA.cu
* Example showing the use of CUFFT for fast 1D-convolution using FFT.
*
* (h_rcvDataIn_short) ==memcpy==>
* (d_signal_padded_float2) ==CUFFTC2C==>
* (d_signal_padded_float2) ==MultByFilter==>
* (d_signal_padded_float2) ==ICUFFTC2C==>
* (d_signal_padded_float2) ==memcpy==>
* (outMatrix)
*
**************************************************************************/

//#define  SHOWALLTIMING 0 //show timing information
//#define  SHOWCOPYTIMING 0 //show memcpy timing information

// The filter size is assumed to be a number smaller than the signal size
#define FILTER_KERNEL_SIZE 12
#define FILTER {-0.0282, 0.0539, 0.1229, -0.2134,  -0.2964, 0.3458, 0.3458, -0.2964, -0.2134, 0.1229, 0.0539, -0.0282}

// includes, system
#include "mex.h"
#include "matrix.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// includes, project
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cufftXt.h>
#include <helper_cuda.h>
#include <helper_functions.h>

// Complex data type
typedef float2 Complex;
static __device__ __host__ inline Complex ComplexScale(Complex, float);
static __device__ __host__ inline Complex ComplexMul(Complex, Complex);
static __global__ void ComplexPointwiseMulAndScale(Complex *, const Complex *, int, float);
static __global__ void floatToComplexPadded(short *, Complex *, int, int, int);
static __global__ void ComplexPaddedToFloat(Complex *, short *, int, int, int);

// Buffers
static Complex *h_padded_filter_kernel=NULL;
static Complex *d_filter_kernel;
short *d_signal;
Complex *d_signal_padded;

// Variables
static cufftHandle plan;
static size_t memSize;
static size_t memSizePadded;
static int new_size;
static int signalSizePadded;
static int Nfft;
static struct cudaPointerAttributes attributes;
static short *outData;  //pointer to output data


// Initialization function
void initialize(int signalSize, int nChs, int numel, short *rcvDataIn);

void cleanup();

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    if (nrhs!=1)
    {
        mexErrMsgIdAndTxt("MyMex:gputest:nrhs", "one inputs required");
    }
    if (!mxIsClass(prhs[0],"int16"))
    {
        mexErrMsgIdAndTxt("MyMex:gputest:prhs", "input must be int16");
    }
    /* General variables */
    mxArray *mxOutArray; //output array (plhs[0])
    short *rcvDataIn;

    const mwSize *dims;
    int signalSize, nChs, numel;
    int const threadsPerBlock = 256;
    int blocksPerGrid;
    /* Setup timing */
    cudaEvent_t start, stop;
    float memcpyH2D_ms, float2complexPadded_ms, FFT_ms, multiply_ms, IFFT_ms, complexPadded2float_ms, memcpyD2H_ms = 0;

    rcvDataIn = mxGetInt16s(prhs[0]); //pointer to incoming data buffer
    numel = mxGetNumberOfElements(prhs[0]);  // get number of elements
    dims = mxGetDimensions(prhs[0]); // data size
    signalSize = dims[0]; // data rows
    nChs = dims[1]; // data cols

    //=== 0. Initialize wall filter and allocate memory the 1st run of the function ===//
    if (d_filter_kernel==NULL)
    {
        initialize(signalSize, nChs, numel, rcvDataIn);
    }

    mxOutArray = mxCreateNumericMatrix(signalSize, nChs, mxINT16_CLASS, mxREAL);  //create output mxArray
    plhs[0] = mxOutArray;
    outData = mxGetInt16s(mxOutArray);

    // Timing objects
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    /*==== 1. Clear padded buffer & Copy RcvData to  Device ====*/
    cudaEventRecord(start,0);
    checkCudaErrors(cudaMemset(d_signal_padded, 0, memSizePadded));  //clear memory on device for padded array

    if( attributes.type == cudaMemoryTypeDevice)
    {
        d_signal = rcvDataIn;  //for GPUDirect: just reassign the pointer
    }
    else
    {
        checkCudaErrors(cudaMemcpyAsync(d_signal, rcvDataIn, memSize, cudaMemcpyHostToDevice,0));
    }
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&memcpyH2D_ms, start, stop);

    /*=== 2. Convert data from float to complex-padded (16bit floating point) & pad data at the same time ===*/
    cudaEventRecord(start,0);
    floatToComplexPadded<<<128,1>>>(d_signal, d_signal_padded, signalSize, signalSizePadded, nChs);
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&float2complexPadded_ms, start, stop);

    /*==== 3. FFT of RcvDataPadded ====*/
    cudaEventRecord(start,0);
    checkCudaErrors(cufftExecC2C(plan, reinterpret_cast<cufftComplex *>(d_signal_padded),
                                 reinterpret_cast<cufftComplex *>(d_signal_padded),
                                 CUFFT_FORWARD));
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&FFT_ms, start, stop);

    /*==== 4.  Multiply RcvDataPadded & FilterKernelPadded ====*/
    cudaEventRecord(start,0);
    blocksPerGrid = (Nfft + threadsPerBlock - 1) / threadsPerBlock;
    ComplexPointwiseMulAndScale<<<blocksPerGrid, threadsPerBlock>>>(d_signal_padded, d_filter_kernel, new_size, (float)1.0f/Nfft);
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&multiply_ms, start, stop);

    /*==== 5. IFFT RcvDataPadded ====*/
    cudaEventRecord(start,0);
    checkCudaErrors(cufftExecC2C(plan, reinterpret_cast<cufftComplex *>(d_signal_padded),
                                 reinterpret_cast<cufftComplex *>(d_signal_padded),
                                 CUFFT_INVERSE));
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&IFFT_ms, start, stop);

    /*==== 6. convert complex-padded to float ====*/
    cudaEventRecord(start,0);
    ComplexPaddedToFloat<<<128,1>>>(d_signal_padded, d_signal, signalSize, signalSizePadded, nChs);
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&complexPadded2float_ms, start, stop);

    /*==== 7. Copy unpadded filtered signal back ====*/
    cudaEventRecord(start,0);
    checkCudaErrors(cudaMemcpyAsync(outData, d_signal, sizeof(short)*numel, cudaMemcpyDeviceToHost,0));
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&memcpyD2H_ms, start, stop);

    #ifdef SHOWALLTIMING
        mexPrintf("MemoryCopy Host-to-Device: %f ms\n"
                  "Float --> Complex-Padded: %f ms\n"
                  "FFT: %f ms\n"
                  "Multiply: %f ms\n"
                  "IFFT: %f ms\n"
                  "ComplexPadded --> Float: %f\n"
                  "MemoryCopy Device-to-Host: %f\n",
                   memcpyH2D_ms, float2complexPadded_ms, FFT_ms, multiply_ms, IFFT_ms, complexPadded2float_ms, memcpyD2H_ms);
    #endif
    #ifdef SHOWCOPYTIMING
        mexPrintf("MemoryCopy Host-to-Device: %f ms\n"
                  "MemoryCopy Device-to-Host: %f ms\n",
                   memcpyH2D_ms, memcpyD2H_ms);
    #endif
}

// Initialize Filters, FFTS, and buffers
void initialize(int signalSize, int nChs, int numel, short* rcvDataIn)
{
    mexPrintf("initialize\n");
    float filter[12] = FILTER;

    Nfft = pow(2, ceil(log2(FILTER_KERNEL_SIZE + (double)signalSize-1))); //calculate optimal size for faster FFTs
    new_size = Nfft*nChs;
    memSize = sizeof(short) * numel;
    memSizePadded = sizeof(Complex) * new_size;
    signalSizePadded = new_size/nChs;

    /* -- Allocate host memory for the signal -- */
    h_padded_filter_kernel = reinterpret_cast<Complex *>(mxCalloc(Nfft*nChs, sizeof(Complex)));

    /* --- Allocate memory on device for filterKernel --- */
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_filter_kernel), memSizePadded));
    checkCudaErrors(cudaMemset(d_filter_kernel, 0, memSizePadded));

    /* --- Allocate a buffer for the incoming RF to be converted to padded float -- */
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_signal), memSize));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_signal_padded), memSizePadded));
    mexMakeMemoryPersistent(h_padded_filter_kernel);

    /* --- Initialize the wall filter ---*/
    for (unsigned int j = 0; j < nChs; ++j)
    {
        for (unsigned int i = 0; i < FILTER_KERNEL_SIZE; ++i)
        {
            h_padded_filter_kernel[i+j*Nfft].x = filter[i];
            h_padded_filter_kernel[i+j*Nfft].y = 0;
        }
    }

    /* --- Copy kernel to device ---*/
    checkCudaErrors(cudaMemcpy(d_filter_kernel, h_padded_filter_kernel, memSizePadded,cudaMemcpyHostToDevice));

    /* --- CUFFT plan simple API ---*/
    checkCudaErrors(cufftPlan1d(&plan, Nfft, CUFFT_C2C, nChs));

    /* --- Copy filter kernel to device ---*/
    checkCudaErrors(cufftExecC2C(plan, reinterpret_cast<cufftComplex *>(d_filter_kernel),
                                       reinterpret_cast<cufftComplex *>(d_filter_kernel),
                                       CUFFT_FORWARD));
    cudaDeviceSynchronize();

    checkCudaErrors(cudaPointerGetAttributes (&attributes, rcvDataIn));

    mexAtExit(cleanup);
}

// Cleanup Function
void cleanup()
{
    if (d_filter_kernel != NULL)
    {
        mexPrintf("Free host memory...");
        mxFree(h_padded_filter_kernel);     //cleanup host memory

        mexPrintf("Destory cuda plans...");
        checkCudaErrors(cufftDestroy(plan));// Destroy CUFFT context

        mexPrintf("Free device memory...");
        cudaFree(d_signal_padded);        // cleanup device memory
        cudaFree(d_filter_kernel);        // cleanup device memory
    }
}

///////////////////////////
// CUDA Kernel Functions //
///////////////////////////


//-- Complex scale --//
static __device__ __host__ inline Complex ComplexScale(Complex a, float s) {
  Complex c;
  c.x = s * a.x;
  c.y = s * a.y;
  return c;
}

//-- Complex multiplication --//
static __device__ __host__ inline Complex ComplexMul(Complex a, Complex b) {
  Complex c;
  c.x = a.x * b.x - a.y * b.y;
  c.y = a.x * b.y + a.y * b.x;
  return c;
}

//-- Complex pointwise multiplication --//
static __global__ void ComplexPointwiseMulAndScale(Complex *a, const Complex *b,
                                                   int size, float scale) {
  const int numThreads = blockDim.x * gridDim.x;
  const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

  for (int i = threadID; i < size; i += numThreads) {
    a[i] = ComplexScale(ComplexMul(a[i], b[i]), scale);
  }
}

//-- Convert float data to complex padded array --//
static __global__ void floatToComplexPadded(short *d_signal, Complex *d_signal_padded,
                                                int signalSize, int signalSizePadded, int nChs)
{
    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
    int stride, paddedStride;
    for (int j = threadID; j < nChs; j += numThreads) //each thread should be aligned to each channel
    {
        paddedStride = j*signalSizePadded;
        stride = j*signalSize;

         for (unsigned int i = 0; i < signalSize; ++i)
         {
             d_signal_padded[i+paddedStride].x = (float)d_signal[i + stride]; //d_signal is padded
             d_signal_padded[i+paddedStride].y = 0; //d_signal is padded
         }
    }
}

//-- Convert complex padded array to float --//
static __global__ void ComplexPaddedToFloat(Complex *d_signal_padded, short *d_signal,
                                                int signalSize, int signalSizePadded, int nChs)
{
    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
    int stride, paddedStride;
    for (int j = threadID; j < nChs; j += numThreads) //each thread should be aligned to each channel
    {
        paddedStride = j*signalSizePadded;
        stride = j*signalSize;
         for (unsigned int i = 0; i < signalSize; ++i)
         {
              d_signal[i + stride] = (short)d_signal_padded[i + paddedStride].x; //d_signal is padded
         }
    }
}
