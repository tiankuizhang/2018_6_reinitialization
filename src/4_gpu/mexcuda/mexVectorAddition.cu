/*==========================================================
 * Filename: mexVectorAddition.cu
 * Description: This function implements element by element
 * vector addition (GPU mex file that contains CUDA code and
 * takes as inputs gpuArray variables)
 * Authors: Ploskas, N., & Samaras, N.
 * Syntax: C = mexVectorAddition(A, b)
 * Input:
 * -- A: a double-precision, floating point GPU array of size N
 * -- B: a double-precision, floating point GPU array of size N
 * Output:
 * -- C: a double-precision, floating point GPU array of size N
 *========================================================*/

#include "mex.h"
#include "gpu/mxGPUArray.h"
/*
 * The kernel code
 */
__global__ void vectorAddition(const double *A,const double *B, double *C, const int N)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i < N){
		C[i] = A[i] + B[i];
	}
}

/*
 * The gateway function
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
	/* declare all variables */
	mxGPUArray const *A;
	mxGPUArray const *B;
	mxGPUArray *C;
	double const *d_A;
	double const *d_B;
	double *d_C;
	int N;

	/* define error messages */
	char const * const errId = "parallel:gpu:mexVectorAddition:InvalidInput";
	char const * const errMsg = "Invalid input to MEX file.";
	
	/* initialize the MathWorks GPU API */
	mxInitGPU();
	
	/* check input data */
	if((nrhs != 2) || !(mxIsGPUArray(prhs[0])) || !(mxIsGPUArray(prhs[1]))) {
		mexErrMsgIdAndTxt(errId, errMsg);
	}

	/* get input arrays */
	A = mxGPUCreateFromMxArray(prhs[0]);
	B = mxGPUCreateFromMxArray(prhs[1]);
	/* verify that A and B are double arrays before extracting the pointer */
	if (mxGPUGetClassID(A) != mxDOUBLE_CLASS || mxGPUGetClassID(B) != mxDOUBLE_CLASS){
		mexErrMsgIdAndTxt(errId, errMsg);
	}

	/* extract a pointer to the input data on the device */
	d_A = (double const *)(mxGPUGetDataReadOnly(A));
	d_B = (double const *)(mxGPUGetDataReadOnly(B));

	/* create a GPUArray to hold the result and get its underlying pointer. */
	C = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(A),
				mxGPUGetDimensions(A),
				mxGPUGetClassID(A),
				mxGPUGetComplexity(A),
				MX_GPU_DO_NOT_INITIALIZE);
	
	d_C = (double *)(mxGPUGetData(C));
	
	/* invoke kernel */
	N = (int)(mxGPUGetNumberOfElements(A));
	int const threadsPerBlock = 256;
	int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
	vectorAddition<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, N);
	
	/* wrap the result up as a MATLAB gpuArray for return */
	plhs[0] = mxGPUCreateMxArrayOnGPU(C);
	
	/*
 	 * the mxGPUArray pointers are host-side structures that refer to
	 * device data. These must be destroyed before leaving the MEX
	 * function.
	 */

	mxGPUDestroyGPUArray(A);
	mxGPUDestroyGPUArray(B);
	mxGPUDestroyGPUArray(C);
}












































