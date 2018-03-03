// reinitializing routine
#include "mex.h"
#include "mexRe.h"
#include "math.h"





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
			mask[ind] = (bool)(lsf[ind]<0) ;
			re_lsf[ind] = lsf[ind];
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
					zpd[up] = dz - zpu[ind]; 
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

	// iteration with runge kutta to update the distance map
	// re_lsf will serve as the current distance map
	double *step = (double *)mxCalloc((mwSize)num_ele, sizeof(double)); // hold time step
	double *intermediate = (double *)mxCalloc((mwSize)num_ele, sizeof(double)); // hold intermediate step

	// 
	for (int i = 0; i < 100; ++i)
	{
		re_step(step, re_lsf, mask, deltat, xpr, xpl, ypf, ypb, zpu, zpd, 
			rows, cols, pges, dx, dy, dz, num_ele);
		for(int j = 0; j < num_ele; ++j)
			intermediate[j] = re_lsf[j] - step[j];
		re_step(step, intermediate, mask, deltat, xpr, xpl, ypf, ypb, zpu, zpd, 
			rows, cols, pges, dx, dy, dz, num_ele);
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

	mxFree(step);
	mxFree(intermediate);

}
