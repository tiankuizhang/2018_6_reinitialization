/*
 * https://research.wmz.ninja/articles/2017/12/using-openmp-in-matlab.html
 * compile on windows: mex -v COMPFLAGS="$COMPFLAGS /openmp" mexAdd.cpp
 * (failed)compile on linux: mex -v CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' mexAdd.cpp 
 * (worked)compile on linux: mex -v mexOpenmpSet.cpp CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"
 */

#include "mex.h"
#include "matrix.h"
#include <omp.h>

/**
 * Call signature: C = mexAdd(A, B)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    // Input validation (omitted)
    // ...
    // ...

    mwSize n1 = mxGetNumberOfElements(prhs[0]);
    mwSize n2 = mxGetNumberOfElements(prhs[1]);
    int nts = (int)mxGetScalar(prhs[2]);
    if (n1 != n2)
    {
        mexErrMsgIdAndTxt("example:add:prhs", "A and B must have the same number of elements.");
    }
    double *A = mxGetPr(prhs[0]);
    double *B = mxGetPr(prhs[1]);
    // Allocate output matrix.
    mxArray *mC = mxCreateDoubleMatrix(n1, 1, mxREAL);
    double *C = mxGetPr(mC);
    
    // Compute the sum in parallel.
    omp_set_num_threads(nts);
    #pragma omp parallel for  default(none) shared(C,A,B,n1) 
    for (int i = 0;i < n1;i++)
    {   
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        //C[i] = A[i] + B[i];
        C[i] = (double)tid;
        mexPrintf("Thread %d : C[%d]=%d \n",tid, i, tid);
    }
    
    // Return the sum.
    plhs[0] = mC;
}