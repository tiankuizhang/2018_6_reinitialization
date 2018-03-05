#include "mex.h"
#include "gpu/mxGPUArray.h"

/*
 * The kernel code
 */

__global__ void 
VectorScalarAddition(const double * A, const double b, double * C, const int N)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	
	if(i < N){
		C[i] = A[i] + b;
	}
}

/*
 * The gateway function
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	const mxGPUArray * A;
	const double * deviceA;
	double b;
	mxGPUArray *C;
	double *deviceC;
	int N;

	char const * const errId = "parallel:gpu:MexVectorScalarAdditionGPU1:InvalidInput";
	char const * const errMsg = "Invalid input to MEX file.";
	
	mxInitGPU();

	if(nrhs !=2 || !(mxIsGPUArray(prhs[0]))){ 
		mexErrMsgIdAndTxt(errId, errMsg);
	}

	A = mxGPUCreateFromMxArray(prhs[0]);
	b = mxGetScalar(prhs[1]);

	if(mxGPUGetClassID(A) != mxDOUBLE_CLASS){
		mexErrMsgIdAndTxt(errId, errMsg);
	}

	deviceA = (double const *)(mxGPUGetDataReadOnly(A));

	C = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(A),
			        mxGPUGetDimensions(A),
				mxGPUGetClassID(A),
				mxGPUGetComplexity(A),
				MX_GPU_DO_NOT_INITIALIZE);
	deviceC = (double *)(mxGPUGetData(C));

	N = (int)(mxGPUGetNumberOfElements(A));
	int const threadsPerBlock = 256;
	int blockPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
	VectorScalarAddition<<<blockPerGrid, threadsPerBlock>>>(deviceA, b, deviceC, N);

	plhs[0] = mxGPUCreateMxArrayOnGPU(C);


	mxGPUDestroyGPUArray(A);
	mxGPUDestroyGPUArray(C);
}
