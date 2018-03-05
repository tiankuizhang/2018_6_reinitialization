__global__ void
matrixMultiplication(const double * A, const double * B, double * C, int N)
{
	int i = blockDim.y * blockIdx.y + threadIdx.y;
	int j = blockDim.x * blockIdx.x + threadIdx.x;
	
	double value = 0;
	
	for(int k = 0; k < N; k++){
		value += A[k * N + j] * B[i * N + k];
	}
	
	C[i * N + j] = value;
}
