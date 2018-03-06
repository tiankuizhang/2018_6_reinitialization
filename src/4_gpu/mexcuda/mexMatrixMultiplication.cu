/*==============================================================================
 *compile with -lcublas option
 *=============================================================================*/
#include "mex.h"
#include "cuda_runtime.h"
#include "cublas_v2.h"

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	
	double * deviceA, * deviceB, * deviceC;
	const double * A, * B;
	double * C;

	int numARows, numACols;
	int numBRows, numBCols;
	int numCRows, numCCols;
	
	char const * const errId = "parallel:gpu:mexMatrixMultiplication:InvalidInput";
	char const * const errMsg = "Invalid input to MEX file.";

	if(nrhs != 2){
		mexErrMsgIdAndTxt(errId, errMsg);
	}

	A = (double *)mxGetData(prhs[0]);
	B = (double *)mxGetData(prhs[1]);


	numARows = (int)mxGetM(prhs[0]);
	numACols = (int)mxGetN(prhs[0]);
	numBRows = (int)mxGetM(prhs[1]);
	numBCols = (int)mxGetN(prhs[1]);

	numCRows = numARows;
	numCCols = numBCols;
	
	plhs[0] = mxCreateNumericMatrix(numCRows, numCCols, mxDOUBLE_CLASS, mxREAL);
	C = (double *)mxGetData(plhs[0]);

	cudaMalloc(&deviceA, sizeof(double) * numARows * numACols);
	cudaMalloc(&deviceB, sizeof(double) * numBRows * numBCols);
	cudaMalloc(&deviceC, sizeof(double) * numCRows * numCCols);

	cublasHandle_t handle;
	cublasCreate(&handle);
	cublasSetMatrix(numARows, 
			numACols,
			sizeof(double),
			A,
			numARows,
			deviceA,
			numARows);
	cublasSetMatrix(numBRows, 
			numBCols,
			sizeof(double),
			B,
			numBRows,
			deviceB,
			numBRows);

	double alpha = 1.0;
	double beta = 0.0;
	cublasDgemm(handle,
		   CUBLAS_OP_N,
		   CUBLAS_OP_N,
		   numARows,
		   numBCols,
		   numACols,
		   &alpha,
		   deviceA,
		   numARows,
		   deviceB,
		   numBRows,
		   &beta,
		   deviceC,
		   numCRows);

	cublasGetMatrix(numCRows,
			numCCols,
			sizeof(double),
			deviceC,
			numCRows,
			C,
			numCRows);

	cublasDestroy(handle);
	cudaFree(deviceA);
	cudaFree(deviceB);
	cudaFree(deviceC);


}

