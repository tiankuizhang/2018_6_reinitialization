__global__ void
vectorAddition(const double * A, const double * B, double * C, const int numElements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if(i < numElements){
		C[i] = A[i] + B[i];
	}
}
