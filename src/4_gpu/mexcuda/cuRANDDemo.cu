/*==========================================================
* Filename: cuRANDDemo.cu
* Description: This function computes the area of the
* definite integral int(x * (x - 2) ^ 6, x = 0..2) (MEX file
* that contains CUDA code and takes as inputs MATLAB arrays)
* Authors: Ploskas, N., & Samaras, N.
* Syntax: area = cuRANDDemo()
* Input: no inputs
* Output:
*
-- area: the area of the definite integral
∗========================================================*/

#include "mex.h"
#include <cuda_runtime.h>
#include <curand.h>
/*
 * The kernel code
 */
__global__ void kernel(int* d_count, double* randomNums)
{
	/* declare all variables */
	int i, tid, xidx, yidx;
	double x, y, z;

	/* Find the overall ID of the thread */
	tid = blockDim.x * blockIdx.x + threadIdx.x;
	i = tid;

	/* Implement the Monte Carlo method */
	xidx = i + i;
	yidx = xidx + 1;

	/*get the random x, y points */
	x = 2 * randomNums[xidx];
	y = 8 * randomNums[yidx];

	/* calculate the value for the selected x point */
	z = x * pow((x - 2), 6);
	/* compare the two values and update the counter */
	if(y <= z){
		d_count[tid] = 1;
	}
	else{
		d_count[tid] = 0;
	}
}

/*
 * The gateway function
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* declare all variables */
	double *B;
	int niter = 1000000;
	double *randomNums;
	double area;
	int threads, blocks;
	int *count, *d_count;
	unsigned int reducedcount = 0;

	/* initialize output array */
	plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
	B = (double *)mxGetData(plhs[0]);

	/* allocate the array for the random numbers */
	cudaMalloc((void**)&randomNums, (2 * niter) * sizeof(double));

	/* use CuRand to generate an array of random numbers on the device */
	curandGenerator_t gen;
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MRG32K3A);
	curandSetPseudoRandomGeneratorSeed(gen, 4294967296ULL^time (NULL));
	curandGenerateUniformDouble(gen, randomNums, 2 * niter);
	curandDestroyGenerator(gen);

	/* set the number of thread blocks and the number of threads per thread block to be launched */
	threads = 1000;
	blocks = 1000;

	/* initialize the array that stores the count values */
	count = (int *)malloc(blocks * threads * sizeof(int));

	/* allocate the array to hold a value (1, 0) whether the point is inside the area (1) or not (0) */
	cudaMalloc((void**)&d_count, blocks * threads * sizeof(int));

	/* launch the kernel */
	kernel<<<blocks, threads>>>(d_count, randomNums);

	/* synchronize the device */
	cudaDeviceSynchronize();

	/* copy the resulting array back */
	cudaMemcpy(count, d_count, blocks * threads * sizeof(int), cudaMemcpyDeviceToHost);
	
	/* Reduce the array into int */
	for(int i = 0; i < niter; i++){
		reducedcount += count[i];
	}

	/* Calculate the area */
	area = ((double)reducedcount / niter) * 2.0 * 8.0;
	
	/* Write area to output */
	*B = area;

	/* destroy the arrays on the device and the host */
	cudaFree(randomNums);
	cudaFree(count);
	free(count);
}
