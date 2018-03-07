/*==============================================================================
 *Filename: cuSPARSEDemo.cu
 * Description: This function implements a sparse matrix-vector
 * multiplication (MEX file that contains CUDA code and takes
 * as inputs MATLAB arrays)
 * Authors: Ploskas, N., & Samaras, N.
 * Syntax: x = cuSPARSEDemo(A, b)
 * Input:
 * -- A: a sparse double-precision, floating point array of size NxN
 * -- b: a double-precision, floating point vector of size Nx1
 * Output:
 * -- x: a double-precision, floating point vector of size Nx1
 *============================================================================*/

#include "mex.h"
#include <cuda_runtime.h>
#include <cusparse.h>

/*
* The gateway function
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	/* declare all variables */
	double *s, *b, *x;
	double *d_s, *d_b, *d_x;
	mwSize *iir, *jjc;
	int *ir, *jc;
	int *d_ir, *d_jc;
	int N, nz;

	/* define error messages */
	char const * const errId = "parallel:gpu:cuSPARSEDemo:InvalidInput";
	char const * const errMsg = "Invalid input to MEX file.";

	/* check input data */
	if(nrhs != 2 || !mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[0]) || mxGetM(prhs[0]) != mxGetN(prhs[0]) || mxGetN(prhs[0]) != mxGetM(prhs[1])){
		mexErrMsgIdAndTxt(errId, errMsg);
	}

	/* initialize output array */
	plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), 1, mxREAL);

	/* get input arrays */
	iir = mxGetIr(prhs[0]); /* Row indexing */
	jjc = mxGetJc(prhs[0]); /* Column count */
	s = mxGetPr(prhs[0]); /* Nonzero elements */
	b = mxGetPr(prhs[1]); /* Rhs vector */
	x = mxGetPr(plhs[0]); /* Output vector */

	/* find arrays dimensions */
	N = mxGetN(prhs[0]);

	/* find the number of nonzero elements */
	nz = (int) mxGetNzmax(prhs[0]);

	/* convert row indexing and column count from mwsize 
	   to int */
	ir = new int[nz];
	jc = new int[N + 1];
	for(int i = 0; i < nz; i++){
		ir[i] = (int)(iir[i]);
	}
	for(int i = 0; i <= N; i++){
	jc[i] = (int)(jjc[i]);
	}

	/* allocate memory on the GPU */
	cudaMalloc((void **)&d_ir, nz * sizeof(int));
	cudaMalloc((void **)&d_jc, (N + 1) * sizeof(int));
	cudaMalloc((void **)&d_s, nz * sizeof(double));
	cudaMalloc((void **)&d_x, N * sizeof(double));
	cudaMalloc((void **)&d_b, N * sizeof(double));
	cudaMemcpy(d_ir, ir, nz * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_jc, jc, (N + 1) * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_s, s, nz * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, N * sizeof(double), cudaMemcpyHostToDevice);

	/* create handle */
	cusparseHandle_t cusparseHandle = 0;
	cusparseCreate(&cusparseHandle);

	/* create and set matrix description and matrix
	   index base */
	cusparseMatDescr_t descr = 0;
	cusparseCreateMatDescr(&descr);
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

	/* perform the sparse matrix-vector multiplication
	   using cuSPARSE */
	double alpha = 1.0;
	double beta = 0.0;
	cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE, N, N, nz, &alpha, descr, d_s, d_jc, d_ir, d_b, &beta, d_x);

	/* copy result */
	cudaMemcpy(x, d_x, N * sizeof(double), cudaMemcpyDeviceToHost);

	/* destroy the handle and the arrays on the device */
	cusparseDestroy(cusparseHandle);
	cudaFree(d_s);
	cudaFree(d_ir);
	cudaFree(d_jc);
	cudaFree(d_b);
	cudaFree(d_x);
}
