#include "mex.h"
#include "cuda_runtime.h"

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
	const double * A;
	double * deviceA;
	double b;
	double *C;
	double *deviceC;
	size_t N;

	char const * const errId = "parallel:gpu:MexVectorScalarAdditionGPU1:InvalidInput";
	char const * const errMsg = "Invalid input to MEX file.";

	if(nrhs !=2 || !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) 
				|| !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
				|| mxGetNumberOfElements(prhs[1]) != 1){
		mexErrMsgIdAndTxt(errId, errMsg);
	}

	A = mxGetPr(prhs[0]);
	b = mxGetScalar(prhs[1]);

	N = mxGetN(prhs[0]);

	plhs[0] = mxCreateDoubleMatrix(1, (mwSize)N, mxREAL);
	C = mxGetPr(plhs[0]);

	cudaMalloc(&deviceA, sizeof(double) * (int)N);
	cudaMalloc(&deviceC, sizeof(double) * (int)N);

	cudaMemcpy(deviceA, A, (int)N * sizeof(double), cudaMemcpyHostToDevice);

	int threadsPerBlock = 256;
	int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
	VectorScalarAddition<<<blocksPerGrid, threadsPerBlock>>>(deviceA, b, deviceC, (int)N);

	cudaMemcpy(C, deviceC, (int)N * sizeof(double), cudaMemcpyDeviceToHost);

	cudaFree(deviceA);
	cudaFree(deviceC);
}

