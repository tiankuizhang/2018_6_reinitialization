/*==============================================================================
 *Filenmae: cuFFTDemo.cu
 * Description: This function implements a two-dimensional
 * discrete Fourier transform of a two-dimensional array
 * using cuFFT library (MEX file that contains CUDA code and
 * takes as inputs MATLAB arrays)
 * Authors: Ploskas, N., & Samaras, N.
 * Syntax: B = cuFFTDemo(A)
 * Input:
 *   -- A: a double-precision, floating point array of size MxN
 * Output:
 *   -- C: a double-precision, floating point array of size 
 *      (M / 2 + 1)/xN
 *============================================================================*/

#include "mex.h"
#include <cuda_runtime.h>
#include <cufft.h>

/*
 * The gateway function
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	/* declare all variables */
	double * deviceA;
	cufftDoubleComplex * deviceB;
	double *A, *B;
	int numARows, numACols;
	int numBRows, numBCols;

	/* define error messages */
	char const * const errId = "parallel:gpu:mexMatrixMultiplication:InvalidInput";
	char const * const errMsg = "Invalid input to MEX file.";
	

	/* check input datea */
	if(nrhs != 1){
		mexErrMsgIdAndTxt(errId, errMsg);
	}

	/* get input array */
	A = (double *)mxGetData(prhs[0]);

	/* find array dimensions */
	numARows = (int)mxGetM(prhs[0]);
	numACols = (int)mxGetN(prhs[0]);

	/* initialize output array */
	numBRows = numARows / 2 + 1;
	numBCols = numACols;

	plhs[0] = mxCreateNumericMatrix(numBRows, numBCols, mxDOUBLE_CLASS, mxREAL);
	B = (double *)mxMalloc(sizeof(cufftDoubleComplex) * numBRows * numBCols);


	/* allocate memory on the GPU */
	cudaMalloc(&deviceA, sizeof(double) * numARows * numACols);
	cudaMalloc(&deviceB, sizeof(cufftDoubleComplex) * numBRows * numBCols);
	cudaMemcpy(deviceA, A, numARows * numACols * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(B, deviceB, numBRows * numBCols * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

	/* create handle and perform the two-dimensional discrete Fourier transform using cuFFT */
	cufftHandle plan;
	cufftPlan2d(&plan, numACols, numARows, CUFFT_D2Z);
	cufftExecD2Z(plan, deviceA, deviceB);

//	cudaMemcpy(B, deviceB, numBRows * numBCols * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

	/* retrieve results */
	double * real = (double *)mxGetPr(plhs[0]);
	double * imag = (double *)mxGetPi(plhs[0]);
	double * complex_ptr = B;
	for(int i = 0; i < numBCols; ++i){
		for(int j = 0; j < numBRows; ++j){
			*real++ = *complex_ptr++;
			*imag++ = *complex_ptr++;
		}
	}

	/* destroy the handle and the arrays on the device */
	cufftDestroy(plan);
	cudaFree(deviceA);
	cudaFree(deviceB);
}
