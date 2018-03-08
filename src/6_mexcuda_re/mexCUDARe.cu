/*******************************************************************************
 * implement reinitialization scheme with CUDA runtime api
 ******************************************************************************/

#include "mex.h"
#include "math.h"
#include "cuda_runtime_api.h"

typedef struct
{
	double sR;
	double sL;
} double_eno_derivative;

__device__ inline
double max2(double x, double y)
{
	return (x<y) ? y : x;
}

__device__ inline
double min2(double x, double y)
{
	return (x<y) ? x : y;
}

__device__ inline
double min_mod(double x, double y)
{
	return (x*y<0) ? 0 : (fabs(x)<fabs(y) ? x : y);
}


__device__ inline
double sign(double x)
{
	return (x>0) ? 1.0 : -1.0;
}


// convert subindex to linear index
// periodic boundary conditions are assumed
__device__ inline
int sub2ind(int row_idx, int col_idx, int pge_idx, int rows, int cols, int pges)
{	
	//int row_idxn = (row_idx + rows) % rows;
	//int col_idxn = (col_idx + cols) % cols;
	//int pge_idxn = (pge_idx + pges) % pges;

	int row_idxn = min2(rows-1, max2(0, row_idx));
	int col_idxn = min2(cols-1, max2(0, col_idx));
	int pge_idxn = min2(pges-1, max2(0, pge_idx));


	int ind = pge_idxn * rows * cols + col_idxn * rows + row_idxn;

	return ind;
}

/****************************************************************************** 
 * calculate distance to the bundary. 
 * if v0*f2<0, return distance from v0 to boundary
 ******************************************************************************/
__device__ inline
double sp(double v1, double v0, double v2, double v3, double ds)
{
	double epsilon=10e-10;
	double p2l = v2 - 2.0 * v0 + v1;
	double p2r = v3 - 2.0 * v2 + v0;
	double p2 = min_mod(p2l, p2r);
	if(p2>epsilon){
		double disc = pow((0.5*p2-v0-v2),2) - 4.*v0*v2;
		double dist =  ds * (0.5 + (v0-v2 - sign(v0-v2)*sqrt(disc)) / p2 );
		return dist;
	}else{
		return ( ds * v0 / (v0 - v2) ); 
	}

}

/****************************************************************************** 
 * calculate Eno derivatives at node v0: [v4,v1,v0,v2,v3]
 ******************************************************************************/
__device__ inline
double_eno_derivative eno_derivative( double v4, double v1, double v0, double v2, double v3, double pr, double pl, double ds)
{
	double p2m;
	double_eno_derivative eno_d;

	double p2 = v1 - 2.0 * v0 + v2;

	double p2r = v0 - 2.0 * v2 + v3;
	p2m = 0.5 * min_mod(p2, p2r) / pow(ds, 2);
	double vr = (pr==ds) ? v2 : 0;
	eno_d.sR = (vr - v0) / pr - pr * p2m;

	double p2l = v0 - 2.0 * v1 + v4;
	p2m = 0.5 * min_mod(p2, p2l) / pow(ds, 2);
	double vl = (pl==ds) ? v1 : 0;
	eno_d.sL = (v0 - vl) / pl + pl * p2m;

	return eno_d;

}

// make corrections to xpr etc
__global__
void boundary_correction(double * xpr, double * xpl, double * ypf, double * ypb, double * zpu, double * zpd, double const * lsf, int num_ele, int rows, int cols, int pges, double dx, double dy, double dz)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	double f2;
	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);
	double f0 = lsf[ind];

	int right = sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
	f2 = lsf[right];
	if(f0*f2<0){
		int left = sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
		int right2 = sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
		xpr[ind] = sp(lsf[left], f0, f2, lsf[right2], dx);
		xpl[right] = dx - xpr[ind];
	}

	int front = sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);
	f2 = lsf[front];
	if(f0*f2<0){
		int back = sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
		int front2 = sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
		ypf[ind] = sp(lsf[back], f0, f2, lsf[front2], dy);
		ypb[front] = dy - ypf[ind];
	}

	int up = sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);
	f2 = lsf[up];
	if(f0*f2<0){
		int down = sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
		int up2 = sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);
		zpu[ind] = sp(lsf[down], f0, f2, lsf[up2], dz);
		zpd[up] = dz - zpu[ind];
	}

}

__global__
void re_step(double * step, double const * lsf, bool const * mask, double const * deltat, double const * xpr, double const * xpl, double const * ypf, double const * ypb, double const * zpu, double const * zpd, int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);

	int right 	= sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
	int right2 	= sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
	int left 	= sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
	int left2 	= sub2ind(row_idx, col_idx-2, pge_idx, rows, cols, pges);

	double_eno_derivative eno_dx = eno_derivative( lsf[left2], lsf[left], lsf[ind], lsf[right], lsf[right2], xpr[ind], xpl[ind], dx);
	double xR = eno_dx.sR;
	double xL = eno_dx.sL;


	int front 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	int front2 	= sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
	int back 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 	= sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);

	double_eno_derivative eno_dy = eno_derivative( lsf[back2], lsf[back], lsf[ind], lsf[front], lsf[front2], ypf[ind], ypb[ind], dy);
	double yF = eno_dy.sR;
	double yB = eno_dy.sL;


	int up 		= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	int up2 	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
	int down 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);

	double_eno_derivative eno_dz = eno_derivative( lsf[down2], lsf[down], lsf[ind], lsf[up], lsf[up2], zpu[ind], zpd[ind], dz);
	double zU = eno_dz.sR;
	double zD = eno_dz.sL;

	if (mask[ind]) {
		step[ind] = ( sqrt(	max2(pow(min2(0,xL),2),pow(max2(0,xR),2)) + 
							max2(pow(min2(0,yB),2),pow(max2(0,yF),2)) + 
							max2(pow(min2(0,zD),2),pow(max2(0,zU),2)) ) - 1)
					* deltat[ind] * (-1.);
	} else{
		step[ind] = ( sqrt(	max2(pow(max2(0,xL),2),pow(min2(0,xR),2)) + 
							max2(pow(max2(0,yB),2),pow(min2(0,yF),2)) + 
							max2(pow(max2(0,zD),2),pow(min2(0,zU),2)) ) - 1)
			* deltat[ind] * (1.);
	}
}

__global__ 
void ini_re(bool * mask, double * deltat, double * xpr, double * xpl, double * ypf, double * ypb, double * zpu, double * zpd, double const * const lsf, int rows, int cols, int pges, double dx, double dy, double dz)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}
	
	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);
	mask[ind] = (lsf[ind]<0);
	deltat[ind] = 0.0;
	xpr[ind] = dx;
	xpl[ind] = dx;
	ypf[ind] = dy;
	ypb[ind] = dy;
	zpu[ind] = dz;
	zpd[ind] = dz;

	
}

__global__
void set_deltat(double * deltat, double const * xpr, double const * xpl, double const * ypf, double const * ypb, double const * zpu, double const * zpd, int rows, int cols, int pges)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}
	
	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);
	double minx = min2(xpr[ind], xpl[ind]);
	double miny = min2(ypf[ind], ypb[ind]);
	double minz = min2(zpu[ind], zpd[ind]);
	deltat[ind] = 0.3 * min2(minx, min2(miny, minz));
}

__global__
void half_step(double * tmp_lsf, double const * step, double const * lsf, int rows, int cols, int pges)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);
	tmp_lsf[ind] = lsf[ind] - step[ind];
}

__global__
void full_step(double * lsf, double const * step, double const * tmp_lsf, int rows, int cols, int pges)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);
	lsf[ind] = (lsf[ind] + tmp_lsf[ind] - step[ind]) / 2;
}

enum in_put {
	level_set_function = 0,
	grid_spacing_x = 1,
	grid_spacing_y = 2,
	grid_spacing_z = 3
};

enum out_put {
	reinitialized_lsf = 0
};

// gateway function
// compile : mexcuda mexCUDARe.cu
// syntax : new_lsf = mexCUDARe(old_lsf, dx, dy, dz);
// takes about 0.05s on a 64x4x64 mesh
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{	
	mxClassID category;

	if(nrhs != 4){
		mexErrMsgIdAndTxt("mexReinitialization:wrong_number_of_inputs", "expecting 4 inputs");
	}

	// assign level set function to LSF
	double *lsf;

	mwSize number_of_dimensions;
	const mwSize *dimension_array;
	size_t num_ele;

	category = mxGetClassID(prhs[level_set_function]);
	number_of_dimensions = mxGetNumberOfDimensions(prhs[level_set_function]);
	dimension_array = mxGetDimensions(prhs[level_set_function]);
	num_ele = mxGetNumberOfElements(prhs[level_set_function]);
	if (category != mxDOUBLE_CLASS || number_of_dimensions != (mwSize)3){
		mexErrMsgIdAndTxt("mexReinitialization:Invalid_Input",
			"Argument %d must be a 3 dimension array of double precision!",
			level_set_function);
	}
	lsf = (double *)mxGetData(prhs[level_set_function]);
	// finish assigning LSF

	// assign size of level set function
	int cols = dimension_array[0];
	int rows = dimension_array[1];
	int pges = dimension_array[2];

	mexPrintf("size of level_set_function: (%d,%d,%d)\n", cols, rows, pges);

	// assign grid spacing array
	double dx = mxGetScalar(prhs[grid_spacing_x]);
	double dy = mxGetScalar(prhs[grid_spacing_y]);
	double dz = mxGetScalar(prhs[grid_spacing_z]);

	mexPrintf("grid_spacing: (%f,%f,%f)\n", dx, dy, dz);

	// create an output array
	double *re_lsf; // pointer to reinitialized level set function
	plhs[reinitialized_lsf] = mxCreateNumericArray(number_of_dimensions, 
		dimension_array, mxDOUBLE_CLASS, mxREAL);
	re_lsf = (double *)mxGetData(plhs[reinitialized_lsf]);
	//  finish creating output array

/*******************************************************************************
 * reinitialization routine begin
 ******************************************************************************/
	int dimx = rows, dimy = 4, dimz = 1;
	dim3 const ThreadBlockSize(dimx, dimy, dimz);
	dim3 const GridSize((rows + ThreadBlockSize.x - 1) / ThreadBlockSize.x, 
		   				(cols + ThreadBlockSize.y - 1) / ThreadBlockSize.y,	
		   				(pges + ThreadBlockSize.z - 1) / ThreadBlockSize.z);	

	// allocate memory for input array
	double * dev_lsf;
	cudaMalloc(&dev_lsf, sizeof(double) * (int)num_ele);
	cudaMemcpy(dev_lsf, lsf, sizeof(double) * (int)num_ele, cudaMemcpyHostToDevice);

	// allocate memory for mask, deltat, xpr etc
	bool * mask;
	double * deltat;
	double * xpr, * xpl, * ypf, * ypb, * zpu, * zpd;
	cudaMalloc(&mask, sizeof(bool) * (int)num_ele);
	cudaMalloc(&deltat, sizeof(double) * (int)num_ele);
	cudaMalloc(&xpr, sizeof(double) * (int)num_ele);
	cudaMalloc(&xpl, sizeof(double) * (int)num_ele);
	cudaMalloc(&ypf, sizeof(double) * (int)num_ele);
	cudaMalloc(&ypb, sizeof(double) * (int)num_ele);
	cudaMalloc(&zpu, sizeof(double) * (int)num_ele);
	cudaMalloc(&zpd, sizeof(double) * (int)num_ele);
	
	ini_re<<<GridSize, ThreadBlockSize>>>(mask, deltat, xpr, xpl, ypf, ypb, zpu, zpd, dev_lsf, rows, cols, pges, dx, dy, dz);

	boundary_correction<<<GridSize, ThreadBlockSize>>>(xpr, xpl, ypf, ypb, zpu, zpd, dev_lsf, num_ele, rows, cols, pges, dx, dy, dz);

	set_deltat<<<GridSize, ThreadBlockSize>>>(deltat, xpr, xpl, ypf, ypb, zpu, zpd, rows, cols, pges);

	double * step, * tmp_lsf;
	cudaMalloc(&step, sizeof(double) * (int)num_ele);
	cudaMalloc(&tmp_lsf, sizeof(double) * (int)num_ele);

	for(int i = 0; i < 100 ; ++i){
		re_step<<<GridSize, ThreadBlockSize>>>(step, dev_lsf, mask, deltat, xpr, xpl, ypf, ypb, zpu, zpd, rows, cols, pges, dx, dy, dz, num_ele);
		half_step<<<GridSize, ThreadBlockSize>>>(tmp_lsf, step, dev_lsf, rows, cols, pges);
		re_step<<<GridSize, ThreadBlockSize>>>(step, tmp_lsf, mask, deltat, xpr, xpl, ypf, ypb, zpu, zpd, rows, cols, pges, dx, dy, dz, num_ele);
		full_step<<<GridSize, ThreadBlockSize>>>(dev_lsf, step, tmp_lsf, rows, cols, pges);
	}

	cudaMemcpy(re_lsf, dev_lsf, sizeof(double) * (int)num_ele, cudaMemcpyDeviceToHost);

	cudaFree(dev_lsf);
	cudaFree(mask);
	cudaFree(deltat);
	cudaFree(xpr);
	cudaFree(xpl);
	cudaFree(ypf);
	cudaFree(ypb);
	cudaFree(zpu);
	cudaFree(zpd);
	cudaFree(step);
	cudaFree(tmp_lsf);
/*******************************************************************************
 * reinitialization routine end
 ******************************************************************************/
}
