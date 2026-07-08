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
* externalDoppler_CPU_MEX.c
* Example showing CDI processing in a mex program
*
*************************************************************************/
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159
#define VELOCITYSCALE (-(253)/(2*PI)) //scale the phase measurement between -pi/pi and color 0 to 255

/* Allocate static Buffers */
static mxComplexDouble *filterMatrix = NULL;
static mxComplexDouble *data;
static mxComplexDouble *dataWF;
static mxComplexDouble *ac;
static double *velocity;
static double *powerSq;
static mxArray *mxCmdPrhs[3];

void Cleanup(void);

void Initialize(int numPixels, int Npri, int nrhs, mxArray const *prhs[]);

void multiplyMatrices(mxComplexDouble *mat1, //matrix 1
                      mxComplexDouble *mat2, //matrix 2
                      mxComplexDouble *res,  //result matrix
                      int nrmat1,   //number of rows for matrix 1
                      int ncmat1,   //number of cols for matrix 1
                      int nrcmat2);  //number of rows & cols for matrix 2 (square)

void lagoneautocorr(mxComplexDouble *IQDataWF,
                    mxComplexDouble *autocorr,
                    mwSize nPixels,
                    mwSize nPri);

void dispMatrix(mxComplexDouble *matrix, int nr, int nc);

void calculatePowerAndVelocity(double *powerSq, double *velocity, mxComplexDouble *ac, mwSize nPixels);

/* Main Function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    // Declare Variables
    int i, j;
    const mwSize *dataSize, *wfDims;
    mwSize dataNumDims, wfNumDims, numPixels, Npri;
    double *IDataIn, *QDataIn, *wfIn;
    const mxArray *threshParams;
    mxArray *imgOut;

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

    // if datasize has changed or was not known before, reinitialize filters
    if (filterMatrix == NULL)
    {
        Initialize((int)numPixels, (int)Npri, nrhs, prhs);
    }

    //////////// Processing ////////////////////
    // copy data in to complex format
    for (i=0; i<(Npri); i++)
    {
        for (j=0; j<(numPixels); j++)
        {
            data[i*numPixels+j].real = IDataIn[i*numPixels+j];
            data[i*numPixels+j].imag = QDataIn[i*numPixels+j];
        }
    }
    /* wallfiltering : multiply data by filterMatrix*/
    multiplyMatrices(data, //matrix 1
                     filterMatrix, //matrix 2
                     dataWF,  //result matrix
                     (int)numPixels,   //number of rows for matrix 1
                     (int)Npri,   //number of cols for matrix 1
                     (int)Npri);  //number of rows & cols for matrix 2 (square)

    /* lag-1 autocorrelation */
    lagoneautocorr(dataWF, ac, numPixels, Npri);

    /* calculate power and velocity*/
    calculatePowerAndVelocity(powerSq, velocity, ac, numPixels);

    /* Use adaptive thresholding */
    threshParams = mexGetVariablePtr("base", "threshParams");
    imgOut = mxCreateNumericMatrix(numPixels, 1, mxDOUBLE_CLASS, mxREAL); // cannot preallocate because of Matlab memory
    mxCmdPrhs[2] = (mxArray *)threshParams;
    //mexCallMATLAB(0, NULL, 0, NULL, "tester");

    mexCallMATLAB(1, &imgOut, 3, mxCmdPrhs, "adaptiveThreshold");
    mxSetDimensions(imgOut, dataSize, dataNumDims-1); // reshape the dimensions for the output matrix
    plhs[0] = imgOut;
}

/*********************************************
*  multiplyMatrices()
*  function to multiply two Complex Matrices
*
**********************************************/

void multiplyMatrices(mxComplexDouble *mat1, //matrix 1
                      mxComplexDouble *mat2, //matrix 2
                      mxComplexDouble *res,  //result matrix
                      int nrmat1,   //number of rows for matrix 1
                      int ncmat1,   //number of cols for matrix 1
                      int nrcmat2)  //number of rows & cols for matrix 2 (square)
{
    int i, j, k;
    for (i = 0; i < ncmat1; i++) //across the row (column index)
    {
        for (j = 0; j < nrmat1; j++) //down the column (row index)
        {
            res[i*nrmat1 + j].real = 0;
            res[i*nrmat1 + j].imag = 0;
            for (k = 0; k < nrcmat2; k++)
            {
                // complex multiply = (x + yi)(u + vi) = (xu – yv) + //real
                //                                       (xv + yu)i  //imag
                res[i*nrmat1+j].real += mat1[k*nrmat1 + j].real*mat2[i*nrcmat2+k].real -
                    mat1[k*nrmat1 + j].imag*mat2[i*nrcmat2+k].imag;
                res[i*nrmat1+j].imag += mat1[k*nrmat1 + j].real*mat2[i*nrcmat2+k].imag +
                    mat1[k*nrmat1 + j].imag*mat2[i*nrcmat2+k].real;
            }
        }
    }

}

/*--- Initialization ---*/
void Initialize(int numPixels, int Npri, int nrhs, mxArray const *prhs[])
{
    int i, ord, wfNumDims;
    double *wfIn;
    const mwSize *wfDims;
    char polystr[80];
    const mxArray *wallFilterWorkspace;

    // Allocate memory
    data = (mxComplexDouble*)(mxCalloc(numPixels*Npri, sizeof(mxComplexDouble)));
    dataWF = (mxComplexDouble*)mxCalloc(numPixels*Npri, sizeof(mxComplexDouble));
    filterMatrix = (mxComplexDouble*)mxCalloc(Npri*Npri, sizeof(mxComplexDouble));
    ac = (mxComplexDouble*)(mxCalloc(numPixels, sizeof(mxComplexDouble)));
    mxCmdPrhs[0] = mxCreateNumericMatrix(numPixels, 1, mxDOUBLE_CLASS, mxREAL); //consider preallocating
    mxCmdPrhs[1] = mxCreateNumericMatrix(numPixels, 1, mxDOUBLE_CLASS, mxREAL); //consider preallocating
    velocity = mxGetDoubles(mxCmdPrhs[0]);
    powerSq = mxGetDoubles(mxCmdPrhs[1]);

    mexMakeMemoryPersistent(data);
    mexMakeMemoryPersistent(dataWF);
    mexMakeMemoryPersistent(filterMatrix);
    mexMakeMemoryPersistent(ac);
//    mexMakeMemoryPersistent(powerSq);
//    mexMakeMemoryPersistent(velocity);
    mexMakeArrayPersistent(mxCmdPrhs[0]);
    mexMakeArrayPersistent(mxCmdPrhs[1]);

    mexAtExit(Cleanup);
    /* -- Wallfilter is an input -- */
    if (nrhs == 3)  //wallfilter has been provided, check size
    {
        wfIn = mxGetDoubles(prhs[2]); // Wallfilter
        wfNumDims = (int)mxGetNumberOfDimensions(prhs[2]);
        wfDims = mxGetDimensions(prhs[2]);
        if (wfDims[0] != Npri)//check size of wallfilter
        {
            mexErrMsgIdAndTxt("MyMex:gputest:wallfilter", "both dimensions of the wallfilter(%i,%i) must equal the PRI(%i) of the ensemble", wfDims[0], wfDims[1], Npri);
        }
        for (i=0; i<(Npri*Npri); i++)
        {
            filterMatrix[i].real = wfIn[i];
            filterMatrix[i].imag = 0;
        }
    }
    else //generate walfilter in Matlab
    {
        // Initialize Wallfilter data using the Verasonics polyfilter.p function
        ord = (int)ceil(Npri/8.0);
        sprintf(polystr, "wallFilterMatrix = polyfilter(%i,%i);", Npri, ord);
        if (!mexEvalString(polystr))
        {
            mexPrintf("Using wall filter created in Matlab\n");
            wallFilterWorkspace = mexGetVariable("base","wallFilterMatrix");
            wfIn = mxGetDoubles(wallFilterWorkspace);
            for (i=0; i<(Npri*Npri); i++)
            {
                filterMatrix[i].real = wfIn[i];
                filterMatrix[i].imag = 0;
            }
        }
    }
}

/*--- Cleanup routine ---*/
void Cleanup(void)
{
    mexPrintf("Clear externalDoppler_CPU_MEX memory\n");
    mxFree(filterMatrix);
    mxFree(data);
    mxFree(dataWF);
    mxFree(ac);
    mxFree(velocity);
    mxFree(powerSq);

}

/*****************************
*  lagoneautocorr()
*  Lag 1 autocorrelation
*
******************************/
void lagoneautocorr(mxComplexDouble *IQDataWF,
                    mxComplexDouble *autocorr,
                    mwSize nPixels,
                    mwSize nPri)
{
    int i, j, n, m;
    double N = (double)(nPri-1);

    for (i=0; i<nPixels; i++)
    {
        autocorr[i].real = 0;
        autocorr[i].imag = 0;
        for (j=0; j<(nPri-1); j++)
        {
            n = j*nPixels;
            m = n + nPixels;
            autocorr[i].real +=  (IQDataWF[i+n].real)*( IQDataWF[i+m].real) +
                (IQDataWF[i+n].imag)*( IQDataWF[i+m].imag);

            autocorr[i].imag +=  (IQDataWF[i+n].real)*( IQDataWF[i+m].imag) -
                (IQDataWF[i+m].real)*( IQDataWF[i+n].imag);
        }
        autocorr[i].real=autocorr[i].real/N; // calculate mean
        autocorr[i].imag=autocorr[i].imag/N; // calculate mean
    }
}

/******************************
 *  calculatePowerAndVelocity()
 *
 ******************************/
void calculatePowerAndVelocity(double *powerSq, double *velocity, mxComplexDouble *ac, mwSize numPixels)
{
    int i;

    for (i=0; i<numPixels; i++)
    {
        powerSq[i] = pow(ac[i].real*ac[i].real + ac[i].imag*ac[i].imag, 0.25);
        velocity[i] = VELOCITYSCALE * atan2(ac[i].imag,ac[i].real);
    }
}