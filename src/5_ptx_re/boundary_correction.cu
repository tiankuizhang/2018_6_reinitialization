/*******************************************************************************
 * make corrections to xpr etc.
 ******************************************************************************/

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
	if(p2>epsilon || p2<-epsilon){
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
