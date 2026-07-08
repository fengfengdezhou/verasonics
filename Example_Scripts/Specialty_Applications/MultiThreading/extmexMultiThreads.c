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
 *
 * extmexMultiThreads.c
 *
 * Template for a multi-threaded external function written as a mex file. This function can be
 * called by runAcq by defining an external processing structure in the user's SetUp script and
 * referencing it in a sequence Event.  The structure for processing IQ data and returning
 * intensity or Doppler data to an ImageBuffer is defined as follows:
 *
 *     Process(1).classname = 'External';
 *     Process(1).method = 'extmexMultiThreads';
 *     Process(1).Parameters = {'srcbuffer','inter',...  % buffer to process, in this case InterBuffer
 *                              'srcbufnum',1,...
 *                              'srcframenum',1,...
 *                              'dstbuffer','image',...  % destination buffer
 *                              'dstbufnum',1,...
 *                              'dstframenum',-1};       % increment frame number for save
 *
 * The function is entered at the mexFunction entry point, where the number of cores is determined and
 * an equal number of threads initialized on the first call.  The thread data structures are then set
 * up and passed to the threads being launched.  The threads do the processing and increment a counter
 * when finished.  The thread that increments the counter to equal the number of cores then is the last
 * thread to finish and signals the mexFunction to continue.  The mexFunction then completes, returning
 * the processing result (if any) to runAcq.
 *
 * You can compile this mex function from the Matlab command line, while in the Vantage directory and
 *  after calling 'activate', using the following commands:
 *
 * >> mex extmexMultiThreads.c -largeArrayDims
 *
 * Before compiling, you will need to define a mex_C_maci64.xml (for Mac) or mexopts.xml (Windows) file that sets the compile parameters.
 *   See provided example files.  Use clang on Mac and VisualStudio on Windows.
 *
 *  Last modified:
 *      10/17/2022
 */

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
#include <sys/sysinfo.h>
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

int getNumCores();

/**** Structure Definitions ****/
// Thread data structure for myThreadFunc. Use this structure to pass parameters to the threads.
//   In this example, the same structure will be passed to all the threads, with the only difference
//   being the threadid.
struct myThreadFuncData {
    int threadid;
    int nrows;
    int ncols;
    // Other variables or pointers.  Any memory used by threads should be pre-allocated and it's
    // pointer passed here.
    //== Edit 3 Start:  Declare the input and output variables ===//
    //short *RFDataIn;
    double *QDataIn;
    double *IDataIn;
    //double *ImgDataIn;

    //short *RFDataOut;
    //double *QDataOut;
    //double *IDataOut;
    double *ImgDataOut;
    //== Edit 3 End ==============================================//

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

//== Edit 1 Start: Define your pixel specific processing function here ===================//
void computeIntensity(double *IDataIn, double *QDataIn, double *ImgDataOut, int pixelIndex)
{
    ImgDataOut[pixelIndex] = sqrt(IDataIn[pixelIndex]*IDataIn[pixelIndex] + QDataIn[pixelIndex]*QDataIn[pixelIndex]);
    return;
}
//== Edit 1 End =========================================================================//

/******************** myThreadFunc Function *****************************/
/*  This is an example thread function used by each thread to compute the magnitude of the IQ pixel
 *  data over a portion of the total pixels, based on the thread ID.
 *  Input: *myThreadFuncData - pointer to a structure of input data with the following attributes
 *                             for this example.
 *    int threadid;
 *    int nrows;
 *    int ncols;
 *    double *QDataIn;
 *    double *IDataIn;
 *    double *ImgDataOut;
 */
void * myThreadFunc(void *myThreadFuncData)
{
    int pixelIndex, id, npixels, nstart, nlimit;
    double x;
    //== Edit 5 Start: Select the needed buffers ===//
    //short *RFDataIn;
    double *IDataIn, *QDataIn;
    //double *ImgDataIn;

    //short *RFDataOut;
    //double *IDataOut, *QDataOUt;
    double *ImgDataOut;
    //== Edit 5 End ================================//

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
        //== Edit 4 Start: non-structure variables for reduced typing ==//
        //RFDataIn = TData->RFDataIn;
        IDataIn = TData->IDataIn;
        QDataIn = TData->QDataIn;
        //ImgDataIn = TData->ImgDataIn;

        //RFDataOut = TData->RFDataOut;
        //IDataOut = TData->IDataOut;
        //QDataOut = TData->QDataOut;
        ImgDataOut = TData->ImgDataOut;
        //== Edit 4 End ================================================//

        /***** Main processing loop for this thread's pixel locations in the InterBuffer. *****/
        // Replace the code below with your own processing routine.
        // - Calculate pixels to process for the threads.
        npixels = TData->nrows * TData->ncols;  // total pixels
        x = (double)npixels/(double)ncores;
        // - calculate nstart and nlimit for this thread.
        nstart = (int)round(x*(double)id);
        nlimit = (int)round(x*(double)(id+1));
        //===============================================================//
        // - Do the I,Q to intensity processing for this threads pixels
        for (pixelIndex=nstart; pixelIndex<nlimit; pixelIndex++) { // for each pixel for this thread

            //== Edit 2 Start: Call your pixel specific function here ==//
            /* This section of code is where the actual pixel specific
            processing occurs.  For this case, i is the pixel index */
            computeIntensity(IDataIn, QDataIn, ImgDataOut, pixelIndex);
            //== Edit 2 End ============================================//
        }
        //===============================================================//
        /***** End of processing for this thread. *****/

        // If threads are being used, increment the finish count and check to see if this thread
        // is the last thread to finish.  If so, signal the main thread and wait to run again.
        #ifndef NOTHREADS
        pthread_mutex_lock(&myThreadFuncCount_mutex);
        pthread_mutex_lock(&myThreadFuncStop_mutex);
        myThreadFuncFinishCount++;
        /* Check the value of count and signal main thread if condition is
           reached.  Note that this occurs while mutex is locked. */
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

    mexPrintf("Cleaning up extmexMultiThreads persistent memory.\n");
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
 *   plhs[0] = ImgDataOut, which is returned at the same dimensions as an ImageBuffer frame.
 */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    int i;
    const mwSize *dims;
    mwSize nRows, nCols;
    //== Edit 6 Select the needed variable ==//
    //short *RFDataIn;
    double *IDataIn, *QDataIn; // note: if the InterBuffer is defined as singles, use "float" datatype
    //double *ImgDataIn;

    //short *RFDataOut;
    //double *IDataOut, *QDataOut; // note: if the InterBuffer is defined as singles, use "float" datatype
    double *ImgDataOut;
    //== Edit 6 End =====================================//
    mxArray *mxa;
    #ifndef NOTHREADS
      int rc;
    #endif

    /*****  Check for proper number of arguments. *****/
    //  NOTE: When using mexErrMsgTxt within an if statement, if true
    //   - mexErrMsgTxt breaks you out of the MEX-file and returns to Matlab.
    if (nrhs != 2) mexErrMsgTxt("Number of inputs should be two - InterBuffer arrays of I & Q data.");
    if (nlhs != 1) mexErrMsgTxt("Number of outputs should be one - array variable for intensity data.");

    // Install cleanup routine to free presistent memory.
    //  - This function will be called by Matlab when executing a 'clear all'.
    mexAtExit(cleanup);


    // Get dimensions from input
    dims = mxGetDimensions(prhs[0]);
    nRows = dims[0];
    nCols = dims[1];



    //== Edit 7 Start: Get pointers from the gateway function to the input and output data ===============//
    // Create the mxArray output array with the size of the InterBuffer and ImageBuffer frame you specified
    //   in your setup script.  If not specified, these buffer frames take the size of PData.
    mxa = mxCreateDoubleMatrix(nRows, nCols, mxREAL);  // memory is allocated here for the output data
    plhs[0] = mxa;
    // uncomment the following if you are returning IQ Data
    //mxb = mxCreateDoubleMatrix(nRows, nCols, mxREAL);  // memory is allocated here for the output data
    //plhs[0] = mxa;

    // Get pointers to the real and imaginary input data.  If the InterBuffer is defined as "single",
    // change the cast for mxGetData.. to (float*).
    //RFDataIn = (short*)mxGetData(prhs[0]);
    QDataIn = (double*)mxGetData(prhs[0]);
    IDataIn = (double*)mxGetData(prhs[1]);
    //ImgDataIn = (double*)mxGetData(prhs[0]);

    //RFDataOut = mxGetPr(mxa);
    //QDataOut = mxGetPr(mxa);
    //IDataOut = mxGetPr(mxb);
    ImgDataOut = mxGetPr(mxa);
    //== Edit 7 End ====================================================================================//

    // If you want to return processed data to the Matlab workspace, you will need to create mxArrays
    // for the data and use mexPutVariable to place the mxArray in the workspace, since any output data
    // will only be returned to the runAcq function, where the extmex function was called.

    /***** Determine number of cores for thread creation. *****/
    #ifdef NOTHREADS
        ncores = 1;
    #else
        ncores = getNumCores();
    #endif

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
        myThreadFuncData[i]->nrows = (int)nRows;
        myThreadFuncData[i]->ncols = (int)nCols;
        //== Edit 8 Start: Data pointers to be passed to each thread ===//
        //myThreadFuncData[i]->RFDataIn = RFDataIn;
        myThreadFuncData[i]->QDataIn = QDataIn;
        myThreadFuncData[i]->IDataIn = IDataIn;
        //myThreadFuncData[i]->ImgDataIn = ImgDataIn;

        //myThreadFuncData[i]->RFDataOut = RFDataOut;
        //myThreadFuncData[i]->QDataOut = QDataOut;
        //myThreadFuncData[i]->IDataOut = IDataOut;
        myThreadFuncData[i]->ImgDataOut = ImgDataOut;
        //== Edit 8 End ================================================//
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

    return;
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
