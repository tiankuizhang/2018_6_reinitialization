#include "mex.h"
#include "omp.h"
#include "math.h"
#include "mexOmpRe.h"

void re_step(double * step, double * re_lsf, bool * mask, double * deltat, 
	double * xpr, double * xpl, double * ypf, double * ypb, double * zpu, double * zpd, 
	int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
{

	/* Adding the following code gives a 4 times speedup!
	 */
	#pragma omp parallel for default(none) shared(step, re_lsf, mask, deltat, \
		xpr, xpl, ypf, ypb, zpu, zpd, rows, cols, pges, dx, dy, dz, num_ele)
	for(int pge_idx = 0; pge_idx < pges; ++pge_idx){
		for(int col_idx = 0; col_idx < cols; ++col_idx){
			for(int row_idx = 0; row_idx < rows; ++row_idx){
				
				int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);	

				int right 	= sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
				int right2 	= sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
				int left 	= sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
				int left2 	= sub2ind(row_idx, col_idx-2, pge_idx, rows, cols, pges);
				

				double_eno_derivative eno_dx = eno_derivative(
					re_lsf[left2], re_lsf[left], re_lsf[ind], re_lsf[right], re_lsf[right2], 
					xpr[ind], xpl[ind], dx);
				double xR = eno_dx.sR;
				double xL = eno_dx.sL;

				
				int front 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
				int front2 	= sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
				int back 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
				int back2 	= sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);

				double_eno_derivative eno_dy = eno_derivative(
					re_lsf[back2], re_lsf[back], re_lsf[ind], re_lsf[front], re_lsf[front2], 
					ypf[ind], ypb[ind], dy);
				double yF = eno_dy.sR;
				double yB = eno_dy.sL;


				int up 		= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
				int up2 	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
				int down 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
				int down2 	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);

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
	}
}