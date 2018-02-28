#include "mex.h"
#include "mexRe.h"
#include "math.h"

void re_step(double * step, double * re_lsf, bool * mask, double * deltat, 
	double * xpr, double * xpl, double * ypf, double * ypb, double * zpu, double * zpd, 
	int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
{

	for(int pge_idx = 0; pge_idx < pges; ++pge_idx){
		for(int col_idx = 0; col_idx < cols; ++col_idx){
			for(int row_idx = 0; row_idx < rows; ++row_idx){
				double f0, f2, p2;
				
				int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);	
				f0 = lsf[ind];					

				int right = sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
				f2 = lsf[right];
				
				
				
			}
		}
	}
}