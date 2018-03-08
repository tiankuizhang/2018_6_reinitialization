/*******************************************************************************
 * more paralle region added. 10 percent faster.
 * source code for reinitialization scheme. combine code in one file for easier
 * organization and precalculate index shift matrices
 * takes 0.23s for a 64x64x64 array. 8 times faster than serial code
 * compile : mex mexOmpReV3.c CFLAGS="$CFLAGS -std=c99 -fopenmp -O3" LDFLAGS="$LDFLAGS -fopenmp"
 ******************************************************************************/

#include "mex.h"
#include "math.h"

typedef struct 
{
	double sR;
	double sL;
} double_eno_derivative;

inline double max2(double x, double y)
{
	return (x<y) ? y : x;
}

inline double min2(double x, double y)
{
	return (x<y) ? x : y;
}

// convert subindex to linear index
// periodic boundary conditions are assumed
inline int sub2ind(int row_idx, int col_idx, int pge_idx, int rows, int cols, int pges)
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

inline double min_mod(double x, double y)
{
	return (x*y<0) ? 0 : (fabs(x)<fabs(y) ? x : y);
}

inline double sign(double x)
{
	return (x>0) ? 1.0 : -1.0;
}

/****************************************************************************** 
 * used in Re.c to calculate distance to the bundary
 * if v0*f2<0, return distance from v0 to boundary
 ******************************************************************************/
inline double sp(double v1, double v0, double v2, double v3, double ds)
{
	double epsilon=10e-10;
	double p2l = v2 - 2.0 * v0 + v1;
	double p2r = v3 - 2.0 * v2 + v0;
	double p2 = min_mod(p2l, p2r);
	if(p2>epsilon){
		double disc = pow((0.5*p2-v0-v2),2) - 4.*v0*v2; 
		double dist = ds * (0.5 + (v0-v2 - sign(v0-v2)*sqrt(disc)) / p2 );
		return dist;
	}else{
		return ( ds * v0 / (v0 - v2) );
	}

}

/****************************************************************************** 
 * used in re_step to calculate reinitialization step
 * calculate Eno derivatives at node v0: [v4,v1,v0,v2,v3]
 ******************************************************************************/
inline double_eno_derivative eno_derivative( double v4, double v1, double v0, double v2, double v3, double pr, double pl, double ds)
{
	double p2m, dsp;
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

void re_step(double * step, double * re_lsf, bool * mask, double * deltat, double * xpr, double * xpl, double * ypf, double * ypb, double * zpu, double * zpd, int * right_ptr, int * left_ptr, int * front_ptr, int * back_ptr, int * up_ptr, int * down_ptr, int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
{
	#pragma omp parallel for default(none) shared(step, re_lsf, mask, deltat, xpr, xpl, ypf, ypb, zpu, zpd, right_ptr, left_ptr, front_ptr, back_ptr, up_ptr, down_ptr, rows, cols, pges, dx, dy, dz, num_ele)
	for(int ind = 0; ind < num_ele; ++ind)
	{
		
		int right 	= right_ptr[ind];
		int right2 	= right_ptr[right_ptr[ind]];
		int left 	= left_ptr[ind];
		int left2 	= left_ptr[left_ptr[ind]];

		double_eno_derivative eno_dx = eno_derivative(
			re_lsf[left2], re_lsf[left], re_lsf[ind], re_lsf[right], re_lsf[right2], 
			xpr[ind], xpl[ind], dx);
		double xR = eno_dx.sR;
		double xL = eno_dx.sL;

		
		int front 	= front_ptr[ind];
		int front2 	= front_ptr[front_ptr[ind]];
		int back 	= back_ptr[ind];
		int back2 	= back_ptr[back_ptr[ind]];

		double_eno_derivative eno_dy = eno_derivative(
			re_lsf[back2], re_lsf[back], re_lsf[ind], re_lsf[front], re_lsf[front2], 
			ypf[ind], ypb[ind], dy);
		double yF = eno_dy.sR;
		double yB = eno_dy.sL;


		int up 		= up_ptr[ind];
		int up2 	= up_ptr[up_ptr[ind]];
		int down 	= down_ptr[ind];
		int down2 	= down_ptr[down_ptr[ind]];

		double_eno_derivative eno_dz = eno_derivative(
			re_lsf[down2], re_lsf[down], re_lsf[ind], re_lsf[up], re_lsf[up2], 
			zpu[ind], zpd[ind], dz);
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
}

void reinitialization(double * re_lsf, double * lsf, int cols, int rows, int pges, 
	double dx, double dy, double dz, int num_ele)
{
	int *right_ptr = 	(int *)mxCalloc((mwSize)num_ele,sizeof(int));
	int *left_ptr = 	(int *)mxCalloc((mwSize)num_ele,sizeof(int));
	int *front_ptr = 	(int *)mxCalloc((mwSize)num_ele,sizeof(int));
	int *back_ptr = 	(int *)mxCalloc((mwSize)num_ele,sizeof(int));
	int *up_ptr = 		(int *)mxCalloc((mwSize)num_ele,sizeof(int));
	int *down_ptr = 	(int *)mxCalloc((mwSize)num_ele,sizeof(int));

	#pragma omp parallel for default(none) shared(right_ptr, left_ptr, front_ptr, back_ptr, up_ptr, down_ptr, pges, cols, rows)
	for(int pge_idx = 0; pge_idx < pges; ++pge_idx){
		for(int col_idx = 0; col_idx < cols; ++col_idx){
			for(int row_idx = 0; row_idx < rows; ++row_idx){
				
				int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);	

				right_ptr[ind] = sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
				left_ptr[ind] = sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);

				front_ptr[ind] = sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
				back_ptr[ind] = sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
				
				up_ptr[ind] = sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
				down_ptr[ind] = sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
			}
		}
	}

	// initialize mask,deltat,xpr etc.
	bool *mask = (bool *)mxCalloc((mwSize)num_ele,sizeof(bool)); // ture for inside, false outside
	double *deltat = (double *)mxCalloc((mwSize)num_ele,sizeof(double)); // time step
	double min_grid = min2(min2(dx,dy),dz);

	double *xpr = (double *)mxCalloc((mwSize)num_ele,sizeof(double)); // distance from boundary : xpr
	double *xpl = (double *)mxCalloc((mwSize)num_ele,sizeof(double));
	double *ypf = (double *)mxCalloc((mwSize)num_ele,sizeof(double));
	double *ypb = (double *)mxCalloc((mwSize)num_ele,sizeof(double));
	double *zpu = (double *)mxCalloc((mwSize)num_ele,sizeof(double));
	double *zpd = (double *)mxCalloc((mwSize)num_ele,sizeof(double));

	#pragma omp parallel for default(none) shared(num_ele, mask, lsf, re_lsf, xpr, xpl, ypf, ypb, zpu, zpd, dx, dy, dz)
	for(int ind = 0; ind < num_ele; ++ind){
			mask[ind] = (bool)(lsf[ind]<0) ;
			re_lsf[ind] = lsf[ind];
			xpr[ind] = dx;
			xpl[ind] = dx;
			ypf[ind] = dy;
			ypb[ind] = dy;
			zpu[ind] = dz;
			zpd[ind] = dz;	
	}

	#pragma omp parallel for default(none) shared(num_ele, lsf, right_ptr, left_ptr, xpr, xpl, dx, front_ptr, back_ptr, ypf, ypb, dy, down_ptr, up_ptr, zpu, zpd, dz)
	for(int ind = 0; ind < num_ele; ++ind){
		
		double f0, f2, p2;
		
		f0 = lsf[ind];					

		int right = right_ptr[ind];
		f2 = lsf[right];
		// if there is a boundary to the right
		if(f0*f2<0){
			int left = left_ptr[ind];
			int right2 = right_ptr[right_ptr[ind]];
			xpr[ind] = sp(lsf[left], f0, f2, lsf[right2], dx);
			xpl[right] = dx - xpr[ind];
		}

		int front = front_ptr[ind];
		f2 = lsf[front];
		// if there is a boundary to the front
		if(f0*f2<0){
			int back = back_ptr[ind];
			int front2 = front_ptr[front_ptr[ind]];
			ypf[ind] = sp(lsf[back], f0, f2, lsf[front2], dy);
			ypb[front] = dy - ypf[ind];		
		}
		
		int up = up_ptr[ind];
		f2 = lsf[up];
		// if there is a boundary to the uppper side
		if(f0*f2<0){
			int down = down_ptr[ind];
			int up2 = up_ptr[up_ptr[ind]];
			zpu[ind] = sp(lsf[down], f0, f2, lsf[up2], dz);
			zpd[up] = dz - zpu[ind]; 
		}
		
	}
	
	#pragma omp parallel for default(none) shared(num_ele, xpr, xpl, ypf, ypb, zpu, zpd, deltat)
	for(int ind = 0; ind < num_ele; ++ind){
		double minx = min2(xpr[ind],xpl[ind]);
		double miny = min2(ypf[ind],ypb[ind]);
		double minz = min2(zpu[ind],zpd[ind]);
		deltat[ind] = 0.3 * min2(minx,min2(miny,minz));	
	}

	// iteration with runge kutta to update the distance map
	// re_lsf will serve as the current distance map
	double *step = (double *)mxCalloc((mwSize)num_ele, sizeof(double)); // hold time step
	double *intermediate = (double *)mxCalloc((mwSize)num_ele, sizeof(double)); // hold intermediate step

	// 
	for (int i = 0; i < 100; ++i)
	{
		re_step(step, re_lsf, mask, deltat, xpr, xpl, ypf, ypb, zpu, zpd, right_ptr, left_ptr, front_ptr, back_ptr, up_ptr, down_ptr, rows, cols, pges, dx, dy, dz, num_ele);
		#pragma omp parallel for default(none) shared(num_ele, intermediate, re_lsf, step)
		for(int j = 0; j < num_ele; ++j)
			intermediate[j] = re_lsf[j] - step[j];
		re_step(step, intermediate, mask, deltat, xpr, xpl, ypf, ypb, zpu, zpd, right_ptr, left_ptr, front_ptr, back_ptr, up_ptr, down_ptr, rows, cols, pges, dx, dy, dz, num_ele);
		#pragma omp parallel for default(none) shared(num_ele, intermediate, re_lsf, step)
		for(int j = 0; j < num_ele; ++j)
			re_lsf[j] = (re_lsf[j] + intermediate[j] - step[j]) / 2;

	}

	mxFree(mask);	
	mxFree(deltat);

	mxFree(xpr);
	mxFree(xpl);
	mxFree(ypf);
	mxFree(ypb);
	mxFree(zpu);
	mxFree(zpd);

	mxFree(right_ptr);
	mxFree(left_ptr);
	mxFree(front_ptr);
	mxFree(back_ptr);
	mxFree(up_ptr);
	mxFree(down_ptr);

	mxFree(step);
	mxFree(intermediate);

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

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{	
	mxClassID category;

	if(nrhs != 4){
		mexErrMsgIdAndTxt("mexReinitialization:wrong_number_of_inputs",
			"expecting 4 inputs");
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

	//	mexPrintf("size of level_set_function: (%d,%d,%d)\n", cols, rows, pges);

	// assign grid spacing array
	double dx = mxGetScalar(prhs[grid_spacing_x]);
	double dy = mxGetScalar(prhs[grid_spacing_y]);
	double dz = mxGetScalar(prhs[grid_spacing_z]);

	// mexPrintf("grid_spacing: (%f,%f,%f)\n", dx, dy, dz);

	// create an output array
	double *re_lsf; // pointer to reinitialized level set function
	plhs[reinitialized_lsf] = mxCreateNumericArray(number_of_dimensions, 
		dimension_array, mxDOUBLE_CLASS, mxREAL);
	re_lsf = (double *)mxGetData(plhs[reinitialized_lsf]);
	//  finish creating output array

	reinitialization(re_lsf, lsf, rows, cols, pges, dx, dy, dz, num_ele);

}
