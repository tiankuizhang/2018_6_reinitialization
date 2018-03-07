/*==========================================================
* Filename: ThrustDemo.cu
* Description: This function calculates the sum of vector
* elements using Thrust (MEX file that contains CUDA code
* and takes as inputs MATLAB arrays)
* Authors: Ploskas, N., & Samaras, N.
* Syntax: b = ThrustDemo(A)
* Input:
*   -- A: a double-precision, floating point vector of size
*	 Nx1 or 1xN
* Output:
*   -- b: the sum of the vector elements
*========================================================*/

#include "mex.h"
#include <thrust/reduce.h>
#include <thrust/device_vector.h>
/*
* The gateway function
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	/* declare all variables */
	double *A, *B;
	int numARows, numACols;
	int N;

	/* define error messages */
	char const * const errId = "parallel:gpu:ThrustDemo:InvalidInput";
	char const * const errMsg = "Invalid input to MEX file.";

	/* check input data */
	if(nrhs != 1 || !mxIsDouble(prhs[0])){
		mexErrMsgIdAndTxt(errId, errMsg);
	}
	
	/* get input array */
	A = (double *)mxGetData(prhs[0]);

	/* find array dimensions */
	numARows = mxGetM(prhs[0]);
	numACols = mxGetN(prhs[0]);
	N = (numARows > numACols) ? numARows : numACols;

	/* initialize output array */
	plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
	B = (double *)mxGetData(plhs[0]);

	/* calculate the sum of the vector using Thurst */
	thrust::device_vector<double> deviceA(A, A + N);
	*B = thrust::reduce(deviceA.begin(), deviceA.end(), 0.0,
		   	thrust::plus<double>());
}
