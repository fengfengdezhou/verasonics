// ----------------------------------------------------------------------
// Copyright 2022 Verasonics, Inc

// Permission is hereby granted, free of charge, to any
// person obtaining a copy of this software and associated
// documentation files (the "Software"), to deal in the
// Software without restriction, including without
// limitation the rights to use, copy, modify, merge,
// publish, distribute, sublicense, and/or sell copies of
// the Software, and to permit persons to whom the Software
// is furnished to do so, subject to the following
// conditions:

// The above copyright notice and this permission notice
// shall be included in all copies or substantial portions
// of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
// ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
// TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
// PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
// SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
// IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
/**************************************************************************
* externalDoppler_CPU_MULTI_MEX.c
* Example showing CDI processing in a mex program with multithreading
*
*************************************************************************/
#define PI 3.14159
#define VELOCITYSCALE ((253)/(2*PI))

// Includes
#include <matrix.h>
#include <math.h>
#include <mat.h>
#include <mex.h>
#include <pthread.h>
#include <matrix.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
#include <windows.h>
#endif
#ifdef __APPLE__
#include <mach/mach.h>  // Includes needed for timer functions. (not currently used)
#include <mach/mach_time.h>
#include <sys/param.h>  // Includes for determining no. of processing cores.
#include <sys/sysctl.h>
#include <unistd.h>
#endif
#ifdef linux
#include <stdint.h>
#include <unistd.h>
#endif
#if defined(_MSC_VER) && defined(__SSE3__)
// Microsoft compiler detected and SSE3 intrinsics requested. SSE3 intrinsics are not used in this
//   example, but could be employed for extra execution speed.
#include    <intrin.h> // For SSE3 intrinsics
#endif
#include <immintrin.h> // For SSE3 (declares _mm_* intrinsics)
#ifdef _WIN32
// Replace unknown round() with equivalent
#define round(value) ((value < 0.0) ? (-(floor(-value + 0.5))) : (floor(value + 0.5)))
#endif
#ifdef _MSC_VER
// Map to equivalent function
#if _MSC_VER<1900  // MSVC 2015 defines C99 compliant snprintf
#define snprintf sprintf_s
#endif
#endif

// Uncomment the #define below if you want to run with a single thread. (This is required if you
//   put in print statements such as mexPrintf() since this print function is not re-entrant.
//#define NOTHREADS

int getNumCores(void);

void getWallFilter(double *wallFilterI, int nPri);

void wallFilter(
        double *dataI,      // data (pixels x pri)
        double *dataQ,      //
        double *filterI,     // filter matrix
        double *dataI_WF,   // filtered data (pixels x pri)
        double *dataQ_WF,   //
        int pixelIndex,
        int nPixels,
        int nPri);

void lagoneautocorr(
        double *dataI,
        double *dataQ,
        double *dataI_AC,
        double *dataQ_AC,
        int pixelIndex,
        int nPixels,
        int nPri);

/**** Structure Definitions ****/
// Thread data structure for myThreadFunc. Use this structure to pass parameters to the threads.
//   In this example, the same structure will be passed to all the threads, with the only difference
//   being the threadid.
struct myThreadFuncData {
    int threadid;
    int nPri;
    int nPixels;
    double pixelsPerCore;
    double *QData;
    double *IData;
};

/**** Global variables ****/
// Most are defined as static, meaning they survive over multiple envocations of the function.
static int ncores, myThreadFuncInit=1;
static volatile int quitThreads=0;
// myThreadFunc Threads globals
static int myThreadFuncThreadsActive=0;
static pthread_t **myThreadFuncThreads;
static struct myThreadFuncData **myThreadFuncData;
int volatile myThreadFuncFinishCount = 0;
static pthread_mutex_t myThreadFuncCount_mutex;
static pthread_cond_t myThreadFuncCount_threshold_cv;
static pthread_mutex_t myThreadFuncStop_mutex;
static pthread_cond_t myThreadFuncStop_cv;

static double* filterI = NULL;
static double* IData_WF;
static double* QData_WF;
static double* IData_AC;
static double* QData_AC;
static double* power;
static double* velocity;
static mxArray *mxCmdPrhs[3];

/******************** myThreadFunc Function *****************************/
/*  This is an example thread function used by each thread to compute the magnitude of the IQ pixel
 *  data over a portion of the total pixels, based on the thread ID.
 *  Input: *myThreadFuncData - pointer to a structure of input data with the following attributes
 *                             for this example.
 *    int threadid;
 *    double pixelsPerCore;
 *    double *QData;
 *    double *IData;
 *    double *ImgData;
 */
void * myThreadFunc(void *myThreadFuncData)
{
    int i, id, nPri, nstart, nlimit, nPixels;
    double pixelsPerCore;
    double *IData, *QData; // pixel x pri
    double *ImgData; // pixel

    struct myThreadFuncData *TData;
    // Processing goes here.  Remember that the addresses for data increments down columns of the
    // Matlab array.  The length of the row and column is not provided with the input, so one needs to either
    // obtain the values from the Matlab workspace (using mexGetVariablePtr), or provide them as #DEFINEs.

    // Initialize structure pointer from input.
    TData = (struct myThreadFuncData *) myThreadFuncData;
    id = TData->threadid;
    // If threads are being used, create an infinite for loop that will execute one pass each time
    // the thread function is called.  Don't use mexPrintf statements if threads are actuve.
    // If NOTHREADS is set, this loop is not needed and taken out during compile.
#ifndef NOTHREADS
    for (;;) {
#endif

        // non-structure variables for reduced typing
        pixelsPerCore = TData->pixelsPerCore;
        nPixels = TData->nPixels;
        nPri = TData->nPri;
        IData = TData->IData;
        QData = TData->QData;


        /***** Main processing loop for this thread's pixel locations in the InterBuffer. *****/
        // Replace the code below with your own processing routine.

        // - calculate nstart and nlimit for this thread.
        nstart = round(pixelsPerCore*(double)id);
        nlimit = round(pixelsPerCore*(double)(id+1));

        // - Do the I,Q to intensity processing for this threads pixels
        for (i=nstart; i<nlimit; i++)  // for each pixel for this thread...
        {
            //1. Wallfilter (multiply ensemble with wallfilter matrix)
            wallFilter(IData, QData, filterI, // filter matrix
                       IData_WF, QData_WF,    // output data
                       i, nPixels, nPri);     // for indexing

            // 2. Lag 1 Autocorrelation
            lagoneautocorr(IData_WF, QData_WF,      // input data
                           IData_AC, QData_AC,     // output data
                           i, nPixels, nPri);      // for indexing

            // 3. Calculate power & velocity
            power[i] = pow((IData_AC[i]*IData_AC[i] + QData_AC[i]*QData_AC[i]),0.25);  //calculate power and add some compression
            velocity[i] = VELOCITYSCALE*(atan2(IData_AC[i],QData_AC[i]));

            // 4.(perform thresholding in Matlab env.)
        }

        /***** End of processing for this thread. *****/

        // If threads are being used, increment the finish count and check to see if this thread
        // is the last thread to finish.  If so, signal the main thread and wait to run again.
#ifndef NOTHREADS
        pthread_mutex_lock(&myThreadFuncCount_mutex);
        pthread_mutex_lock(&myThreadFuncStop_mutex);
        myThreadFuncFinishCount++;
        /* Check the value of count and signal main thread if condition is
         * reached.  Note that this occurs while mutex is locked. */
        if (myThreadFuncFinishCount == ncores)
            pthread_cond_signal(&myThreadFuncCount_threshold_cv);
        pthread_mutex_unlock(&myThreadFuncCount_mutex);

        /* Wait for signal from main thread to run again. */
        pthread_cond_wait(&myThreadFuncStop_cv, &myThreadFuncStop_mutex);
        pthread_mutex_unlock(&myThreadFuncStop_mutex);
        if (quitThreads==1) break;
    }
#endif
    return NULL;
}

/******************** Cleanup Function  ********************/
/* This function is called when the MEX file is cleared by Matlab (as with a 'clear all'.
 * It is used to free memory allocated for thread parameters and persistent memory. */

void cleanup(void) {
    int i;

    mexPrintf("Cleaning up externalDoppler_CPU_MEX_MULTI persistent memory.\n");
    /* Clean up active threads for xmitField processing. */
    quitThreads = 1;
    if (myThreadFuncThreadsActive != 0) {
        /* Release the myThreadFunc threads so they can exit. */
        pthread_mutex_lock(&myThreadFuncStop_mutex);
        pthread_cond_broadcast(&myThreadFuncStop_cv);
        pthread_mutex_unlock(&myThreadFuncStop_mutex);
        for (i=0; i<ncores; i++) {
            pthread_join(*myThreadFuncThreads[i], NULL);
        }
        pthread_mutex_destroy(&myThreadFuncCount_mutex);
        pthread_cond_destroy(&myThreadFuncCount_threshold_cv);
        pthread_mutex_destroy(&myThreadFuncStop_mutex);
        pthread_cond_destroy(&myThreadFuncStop_cv);
    }
    if (myThreadFuncInit == 0) {
        // Free memory for thread structures.
        mxFree(myThreadFuncThreads[0]);
        mxFree(myThreadFuncThreads);
        mxFree(myThreadFuncData[0]);
        mxFree(myThreadFuncData);
    }

    free(IData_WF);
    free(QData_WF);
    free(IData_AC);
    free(QData_AC);
    mxFree(velocity);
    mxFree(power);
}


/********** mexFunction - The MATLAB external function call entry point. *******
 *
 * When MATLAB calls a C function, it defines the following:
 *   [a,b,c,...] = fun(d,e,f,...)
 * The items on the left (a,b,c) are referred to as "left hand side" arguments
 * and are the returned items from this function.  Likewise, (d,e,f) are
 * "right hand side" arguments and refer to inputs to this function.
 *
 * Input/Output
 *   nlhs - "Number of Left Hand Side" parameters.
 *   plhs - nlhs sized array of pointers that point to mxArray outputs
 *   nrhs - "Number of Right Hand Side" parameters.
 *   prhs - nrhs sized array of pointers that point to mxArray inputs
 *
 * For this function, the following inputs are required:
 *   nlhs = 1;
 *   prhs[0] - InterBuffer complex mexArray pointer.
 *
 * Output returned:
 *   nlhs = 1;
 *   plhs[0] = myImgData, which is returned at the same dimensions as an ImageBuffer frame.
 */

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    int i;
    int nRows, nCols, nPixels, nPri, dataNumDims;
    const mwSize *wfDims;
    const mwSize *dataDims;
    double pixelsPerCore; //note: this may be a non-integer
    double *myIData, *myQData; // note: if the InterBuffer is defined as singles, use "float" datatype
    mxArray *imgOut;
    const mxArray *threshParams;

#ifndef NOTHREADS
    int rc;
#endif

    /*****  Check for proper number of arguments. *****/
    //  NOTE: When using mexErrMsgTxt within an if statement, if true
    //   - mexErrMsgTxt breaks you out of the MEX-file and returns to Matlab.
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
    dataDims = mxGetDimensions(prhs[0]);
    nRows = dataDims[0];
    nCols = dataDims[1];
    if (dataNumDims == 4)
    {
        nPixels = nRows*nCols*dataDims[2];
        nPri = dataDims[3];
    }
    else if (dataNumDims == 3)
    {
        nPixels = nRows*nCols;
        nPri = dataDims[2];
    }
    else
    {
        mexErrMsgTxt("This function can only process data with 3 or 4 dimensions.");
    }
    if (filterI == NULL)
    {
        mexPrintf("Initialize memory\n");
        filterI = calloc(nPri*nPri, sizeof(double));
        IData_WF = calloc(nPixels*nPri, sizeof(double));   //temp  need to free
        QData_WF = calloc(nPixels*nPri, sizeof(double));   //temp
        IData_AC = calloc(nPixels, sizeof(double));   //temp  need to free
        QData_AC = calloc(nPixels, sizeof(double));   //temp
        power = calloc(nPixels, sizeof(double));   //temp  need to free
        mxCmdPrhs[0] = mxCreateNumericMatrix(nPixels, 1, mxDOUBLE_CLASS, mxREAL);
        mxCmdPrhs[1] = mxCreateNumericMatrix(nPixels, 1, mxDOUBLE_CLASS, mxREAL);
        velocity = mxGetDoubles(mxCmdPrhs[0]);
        power = mxGetDoubles(mxCmdPrhs[1]);

        mexMakeMemoryPersistent(filterI);
        mexMakeMemoryPersistent(IData_WF);
        mexMakeMemoryPersistent(QData_WF);
        mexMakeMemoryPersistent(IData_AC);
        mexMakeMemoryPersistent(QData_AC);
        mexMakeMemoryPersistent(power);
        mexMakeMemoryPersistent(velocity);
        mexMakeArrayPersistent(mxCmdPrhs[0]);
        mexMakeArrayPersistent(mxCmdPrhs[1]);

        if (nrhs == 3)  //wallfilter has been provided, check size
        {
            filterI = mxGetDoubles(prhs[2]); // Wallfilter
            wfDims = mxGetDimensions(prhs[2]);
            if (wfDims[0] != nPri) // check size of wallfilter
            {
                mexErrMsgIdAndTxt("MyMex:gputest:wallfilter", "both dimensions of the wallfilter(%i,%i) must equal the PRI(%i) of the ensemble", wfDims[0], wfDims[1], nPri);
            }
        }
        else //generate wall filter
        {
            getWallFilter(filterI, nPri);
        }
    }

    // Install cleanup routine to free presistent memory.
    //  - This function will be called by Matlab when executing a 'clear all'.
    mexAtExit(cleanup);

    // Get pointers to the real and imaginary input data.  If the InterBuffer is defined as "single",
    // change the cast for mxGetData.. to (float*).
    myQData = (double*)mxGetData(prhs[0]);
    myIData = (double*)mxGetData(prhs[1]);

    /***** Determine number of cores for thread creation. *****/
#ifdef NOTHREADS
    ncores = 1;
#else
    ncores = getNumCores();
#endif
    pixelsPerCore = (double)nPixels/(double)ncores;
    /***** Initialization. *****/

    // If executing for the first time, perform some initialization functions.
    if (myThreadFuncInit == 1) {  // myThreadFuncInit is static variable and persists over multiple calls
        // Allocate memory for thread structures and myThreadFuncData structures
        if ((myThreadFuncThreads = mxCalloc(ncores,sizeof(int *)))) mexMakeMemoryPersistent(myThreadFuncThreads);
        else mexErrMsgTxt("extmexMultiThreads: out of memory error.\n");
        if ((myThreadFuncThreads[0] = mxCalloc(ncores, sizeof(pthread_t)))) mexMakeMemoryPersistent(myThreadFuncThreads[0]);
        else mexErrMsgTxt("extmexMultiThreads: out of memory error.\n");
        if ((myThreadFuncData = mxCalloc(ncores,sizeof(int *)))) mexMakeMemoryPersistent(myThreadFuncData);
        else mexErrMsgTxt("extmexMultiThreads: out of memory error.\n");
        if ((myThreadFuncData[0] = mxCalloc(ncores, sizeof(struct myThreadFuncData)))) mexMakeMemoryPersistent(myThreadFuncData[0]);
        else mexErrMsgTxt("extmexMultiThreads: out of memory error.\n");
        for (i=0; i<ncores; i++) {
            myThreadFuncThreads[i] = myThreadFuncThreads[0] + i;
            myThreadFuncData[i] = myThreadFuncData[0] + i;
        }
        myThreadFuncInit = 0;
    }

    /***** Fill in myThreadFuncData parameters and launch the threads. *****/

    for (i=0; i<ncores; i++) {
        /* Put all the parameters in each ThreadFuncData structure. */
        myThreadFuncData[i]->threadid = i;
        myThreadFuncData[i]->nPixels = nPixels;
        myThreadFuncData[i]->pixelsPerCore = pixelsPerCore;
        myThreadFuncData[i]->nPri = nPri;
        myThreadFuncData[i]->QData = myQData;
        myThreadFuncData[i]->IData = myIData;
    }
#ifdef NOTHREADS
    myThreadFunc(myThreadFuncData[0]);
#else
    if (myThreadFuncThreadsActive == 0) {
        /* Initialize mutex and condition variable objects */
        pthread_mutex_init(&myThreadFuncCount_mutex, NULL);
        pthread_cond_init (&myThreadFuncCount_threshold_cv, NULL);
        pthread_mutex_init(&myThreadFuncStop_mutex, NULL);
        pthread_cond_init (&myThreadFuncStop_cv, NULL);
        myThreadFuncFinishCount = 0;
        /* Launch the threads. */
        for (i=0; i<ncores; i++) {
            rc = pthread_create(myThreadFuncThreads[i], NULL, myThreadFunc, (void *)myThreadFuncData[i]);
            if (rc) mexErrMsgTxt("Problem launching myThreadFunc threads.\n");
        }
        //mexPrintf("Launched the threads.\n");
        myThreadFuncThreadsActive = 1;
    }
    else {
        /* Release the threads to run again. */
        myThreadFuncFinishCount = 0;
        pthread_mutex_lock(&myThreadFuncStop_mutex);
        pthread_cond_broadcast(&myThreadFuncStop_cv);
        pthread_mutex_unlock(&myThreadFuncStop_mutex);
    }

    pthread_mutex_lock(&myThreadFuncCount_mutex);
    //mexPrintf("xmitFieldFinishCount = %d\n",xmitFieldFinishCount);
    if (myThreadFuncFinishCount<ncores) {
        pthread_cond_wait(&myThreadFuncCount_threshold_cv, &myThreadFuncCount_mutex);
        //mexPrintf("myThreadFuncFinishCount = %d: thread Condition signal received.\n",myThreadFuncFinishCount);
    }
    pthread_mutex_unlock(&myThreadFuncCount_mutex);
#endif

    // The intensity data computed by the threads was stored in the memory of the output array phls[0]
    // so there is no need for further processing by the mexFunction.
    imgOut = mxCreateNumericMatrix(nPixels, 1, mxDOUBLE_CLASS, mxREAL); // cannot preallocate because of Matlab memory
    threshParams = mexGetVariablePtr("base", "threshParams");
    if (threshParams == NULL) mexErrMsgTxt("threshParams must be defined in base workspace");
    mxCmdPrhs[2] = (mxArray *)threshParams;
    mexCallMATLAB(1, &imgOut, 3, mxCmdPrhs, "adaptiveThreshold");
    mxSetDimensions(imgOut, dataDims, dataNumDims-1); // reshape the dimensions for the output matrix
    plhs[0] = imgOut;
    return;
}

/*****************************
 * getWallFilter()
 *
 *****************************/
void getWallFilter(double *filterI, int nPri)
{
    int i, ord;
    char polystr[80];
    mxArray *wallFilterWorkspace;
    mxArray *polyfilterInputs[2];
    mxDouble *filterPtr;

    // Initialize wallfilter using the Verasonics polyfilter.p function
    ord = (int)ceil(nPri/8.0);
    polyfilterInputs[0] = mxCreateDoubleScalar(nPri);
    polyfilterInputs[1] = mxCreateDoubleScalar(ord);

    if (!mexCallMATLAB(1, &wallFilterWorkspace, 2, polyfilterInputs, "polyfilter"))
    {
    //mexPutVariable("base", "TESTVAR", wallFilterWorkspace); //for debugging
    filterPtr = mxGetDoubles(wallFilterWorkspace);

    // Copy filter from base workspace to static memory
    memcpy(filterI, filterPtr, nPri*nPri*sizeof(double));
    }
    else
    {
        mexErrMsgTxt("polyfilter error");
    }
}

/*****************************
 *  lagoneautocorr()
 *  Lag 1 autocorrelation
 *
 ******************************/
void lagoneautocorr(
        double *dataI,
        double *dataQ,
        double *dataI_AC,
        double *dataQ_AC,
        int pixelIndex,
        int nPixels,
        int nPri)
{
    int j, n, m;
    double N = (double)(nPri-1);

    dataI_AC[pixelIndex] = 0;
    dataQ_AC[pixelIndex] = 0;
    for (j=0; j<(nPri-1); j++)
    {
        n = j*nPixels;
        m = n + nPixels;
        dataI_AC[pixelIndex] += (dataI[pixelIndex + n])*(dataI[pixelIndex+m]) +
                (dataQ[pixelIndex + n])*(dataQ[pixelIndex+m]);

        dataQ_AC[pixelIndex] += (dataI[pixelIndex+n])*(dataQ[pixelIndex+m]) -
                (dataI[pixelIndex+m])*(dataQ[pixelIndex+n]);
    }
    dataI_AC[pixelIndex] = dataI_AC[pixelIndex]/N; //calculate mean
    dataQ_AC[pixelIndex] = dataQ_AC[pixelIndex]/N; //calculate mean
}


/*********************************************
 *  wallFilter()
 *  function to multiply two Complex Matrices
 *
 *********************************************/
void wallFilter(
        double *dataI,      // data (pixels x pri)
        double *dataQ,      //
        double *filterI,    // filter matrix
        double *dataI_WF,   // filtered data (pixels x pri)
        double *dataQ_WF,   //
        int pixelIndex,
        int nPixels,
        int nPri)
{
    int j, k;
    int n;

    for (j = 0; j< nPri; j++)
    {
        dataI_WF[pixelIndex + j*nPixels] = 0.0;
        dataQ_WF[pixelIndex + j*nPixels] = 0.0;

        for (k = 0; k < nPri; k++)
        {
            // complex multiply = (x + yi)(u + vi) = (xu – yv) + // real
            //                                       (xv + yu)i  // imag
            //dataI is "2D" pixel x pri

            dataI_WF[pixelIndex + j*nPixels] += dataI[pixelIndex  + k*nPixels] * filterI[j*nPri + k];  // imag. part of WF is zero
            dataQ_WF[pixelIndex + j*nPixels] += dataQ[pixelIndex  + k*nPixels] * filterI[j*nPri + k];  // imag. part of WF is zero
            //mexPrintf("= %f\n", dataI_WF[pixelIndex + j*nPixels]);
        }

    }
}

// C utility function to get the number of logical cores.  Note that this is not the same as
//   the number of physical cores if the processors are hyperthreading capable.
int getNumCores() {
    int nCores;
#ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    nCores = sysinfo.dwNumberOfProcessors;
#elif __APPLE__
    int nm[2];
    size_t len = 4;
    uint32_t count;

    nm[0] = CTL_HW; nm[1] = HW_AVAILCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);

    if(count < 1) {
        nm[1] = HW_NCPU;
        sysctl(nm, 2, &count, &len, NULL, 0);
        if(count < 1) { count = 1; }
    }
    nCores = count;
#else
    nCores = sysconf(_SC_NPROCESSORS_ONLN);
#endif
    nCores = nCores/2; //if hyperthreading
    return nCores;
}
