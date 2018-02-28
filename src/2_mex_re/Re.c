// reinitializing routine
#include "mex.h"
#include "mexRe.h"
#include "math.h"


// convert subindex to linear index
// periodic boundary conditions are assumed
inline int sub2ind(int row_idx, int col_idx, int pge_idx, int rows, int cols, int pges)
{	
	int row_idxn = (row_idx + rows) % rows;
	int col_idxn = (col_idx + cols) % cols;
	int pge_idxn = (pge_idx + pges) % pges;

	int ind = pge_idxn * rows * cols + col_idxn * rows + row_idxn;

	return ind;
}

inline double max2(double x, double y)
{
	return (x<y) ? y : x;
}

inline double min2(double x, double y)
{
	return (x<y) ? x : y;
}

inline double min_mod(double x, double y)
{
	return (x*y<0) ? 0 : (fabs(x)<fabs(y) ? x : y);
}

inline double sign(double x)
{
	return (x<0) ? -1.0 : 1.0;
}

// if f0*f2<0, return distance from f0 to boundary
inline double sp(double v1, double v0, double v2, double v3, double ds)
{
	double epsilon=10e-10;
	double p2l = v2 - 2.0 * v0 + v1;
	double p2r = v3 - 2.0 * v2 + v0;
	double p2 = min_mod(p2l, p2r);
	if(p2>epsilon){
		return dist(Discriminant(p2,v0,v2),ds,p2,v0,v2);
	}else{
		return dist_turn(ds,v0,v2);
	}

}

inline double Discriminant(double p2, double v0, double v2)
{
	return ( pow((0.5*p2-v0-v2),2) - 4.*v0*v2 );
}

inline double dist(double disc, double ds,double p2,double v0, double v2)
{
	return ( ds * (0.5 + (v0-v2 - sign(v0-v2)*sqrt(disc)) / p2 ) );
}

inline double dist_turn(double ds, double v0, double v2)
{
	return ( ds * v0 / (v0 - v2) );	
}
// above are some scaffold functions to calculate sp



void reinitialization(double * re_lsf, double * lsf, int cols, int rows, int pges, 
	double dx, double dy, double dz, int num_ele)
{
	mexPrintf("reinitialization routine started \n");

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

	for(int ind = 0; ind < num_ele; ++ind){
			mask[ind] = lsf[ind]<0 ;
			xpr[ind] = dx;
			xpl[ind] = dx;
			ypf[ind] = dy;
			ypb[ind] = dy;
			zpu[ind] = dz;
			zpd[ind] = dz;	
	}

	for(int pge_idx = 0; pge_idx < pges; ++pge_idx){
		for(int col_idx = 0; col_idx < cols; ++col_idx){
			for(int row_idx = 0; row_idx < rows; ++row_idx){
				double f0, f2, p2;
				
				int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);	
				f0 = lsf[ind];					

				int right = sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
				f2 = lsf[right];
				// if there is a boundary to the right
				if(f0*f2<0){
					int left = sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
					int right2 = sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
					xpr[ind] = sp(lsf[left], f0, f2, lsf[right2], dx);
					xpl[right] = dx - xpr[ind];
				}

				int front = sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
				f2 = lsf[front];
				// if there is a boundary to the front
				if(f0*f2<0){
					int back = sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
					int front2 = sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
					ypf[ind] = sp(lsf[back], f0, f2, lsf[front2], dy);
					ypb[front] = dy - ypf[ind];		
				}
				
				int up = sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
				f2 = lsf[up];
				// if there is a boundary to the uppper side
				if(f0*f2<0){
					int down = sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
					int up2 = sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);		
					zpu[ind] = sp(lsf[down], f0, f2, lsf[up2], dz);
					zpd[ind] = dz - zpu[ind];
				}
				
			}
		}
	}

	for(int ind = 0; ind < num_ele; ++ind){
		double minx = min2(xpr[ind],xpl[ind]);
		double miny = min2(ypf[ind],ypb[ind]);
		double minz = min2(zpu[ind],zpd[ind]);
		deltat[ind] = 0.3 * min2(minx,min2(miny,minz));	
	}

	// iteraction with runge kutta to update the distance map
	// re_lsf will serve as the current distance map
	double *step = (double *)mxCalloc((mwSize)num_ele, sizeof(double)); // hold time step
	double *intermediate = (double *)mxCalloc((mwSize)num_ele, sizeof(double)); // hold intermediate step







	mxFree(mask);	
	mxFree(deltat);
	mxFree(xpr);
	mxFree(xpl);
	mxFree(ypf);
	mxFree(ypb);
	mxFree(zpu);
	mxFree(zpd);
	mxFree(step);
	mxFree(intermediate);

}
