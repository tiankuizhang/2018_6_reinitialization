/*==============================================================================
 *Filename: cuSOLVERDemo.cu
 *Description: This function solves a system of linear
 *equations Ax = b for x (MEX file that contains CUDA code
 *and takes as inputs MATLAB arrays)
 *Authors: Ploskas, N., & Samaras, N.
 *Syntax: x = cuSOLVERDemo(A, b)
 *Input:
 *  -- A: a double-precision, floating point array of size NxN
 *  -- b: a double-precision, floating point vector of size Nx1
 * Output:
 *	-- x: a double-precision, floating point vector of size Nx1
 *============================================================================*/

#include "mex.h"
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

/*
 * The gateway function
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* declare all variables */
	double *deviceA, *deviceB, *deviceX;
	double *A, *b, *x;
	int N;
	int lwork = 0;
	double *d_work = NULL;
	int *devInfo = NULL;
	const double one = 1;

	/* define error messages */
	char const * const errId = "parallel:gpu:cuSOLVERDemo:InvalidInput";
	char const * const errMsg = "Invalid input to MEX file.";

	/* check input data */
	if(nrhs != 2 || !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || mxGetM(prhs[0]) != mxGetN(prhs[0]) || mxGetM(prhs[0]) != mxGetM(prhs[1])){
		mexErrMsgIdAndTxt(errId, errMsg);
	}

	/* get input arrays */
	A = (double *)mxGetData(prhs[0]);
	b = (double *)mxGetData(prhs[1]);
	
	/* find arrays dimensions */
	N = (int)mxGetN(prhs[0]);
	
	/* initialize output array */
	plhs[0] = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);
	x = (double *)mxGetData(plhs[0]);

	/* allocate memory on the GPU */
	cudaMalloc(&deviceA, sizeof(double) * N * N);
	cudaMalloc(&deviceB, sizeof(double) * N);
	cudaMalloc(&deviceX, sizeof(double) * N);
	cudaMemcpy(deviceA, A, sizeof(double) * N * N, cudaMemcpyHostToDevice);
	cudaMemcpy(deviceB, b, sizeof(double) * N, cudaMemcpyHostToDevice);
	
	/* create handle */
	cusolverDnHandle_t cudenseH = NULL;
	cusolverDnCreate(&cudenseH);

	/* query working space of geqrf */
	cusolverDnDgeqrf_bufferSize(cudenseH, N, N, deviceA, N, &lwork);
	cudaMalloc((void**)&d_work, sizeof(double) * lwork);
	cudaMalloc((void **)&devInfo, sizeof(int));

	/* compute QR factorization */ 
	cusolverDnDgeqrf(cudenseH, N, N, deviceA, N, deviceX, d_work, lwork, devInfo);

	/* synchronize the device */
	cudaDeviceSynchronize();

	/* compute Q^T*b */
	cusolverDnDormqr(cudenseH, CUBLAS_SIDE_LEFT, CUBLAS_OP_T, N, 1, N, deviceA, N, deviceX, deviceB, N, d_work, lwork, devInfo);

	/* synchronize the device */
	cudaDeviceSynchronize();

	/* create handle */
	cublasHandle_t cublasH = NULL;
	cublasCreate(&cublasH);

	/* compute x = R \ Q^T*b */
	cublasDtrsm(cublasH, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, N, 1, &one, deviceA, N, deviceB, N);

	/* synchronize the device */
	cudaDeviceSynchronize();

	/* copy the result */
	cudaMemcpy(x, deviceB, sizeof(double) * N, cudaMemcpyDeviceToHost);

	/* destroy the handles and the arrays on the device */
	cublasDestroy(cublasH);
	cusolverDnDestroy(cudenseH);
	cudaFree(deviceA);
	cudaFree(deviceB);
	cudaFree(deviceX);
}
