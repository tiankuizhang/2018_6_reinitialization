#include "mex.h"
#include "mexRe.h"
#include "math.h"


// calculate Eno derivatives at node v0: [v4,v1,v0,v2,v3]
inline double eno_derivative(double v4, double v1, double v0, double v2, double v3, 
	double pr, double pl, double ds, double * sR, double * sL)
{
	double p2m, dsp;

	double p2 = v1 - 2.0 * v0 + v2;

	double p2r = v0 - 2.0 * v2 + v3;
	p2m = 0.5 * min_mod(p2, p2r) / pow(ds, 2);
	dsp = (pr<ds) ? pr : ds;
	double vr = (pr<ds) ? 0 : v2;
	(*sR) = (vr - v0) / dsp - dsp * p2m;

	double p2l = v0 - 2.0 * v1 + v4;
	p2m = 0.5 * min_mod(p2, p2l) / pow(ds, 2);
	dsp = (pl<ds) ? pl : ds;
	double vl = (pl<ds) ? 0 : v1;
	(*sL) = (v0 - vl) / dsp + dsp * p2m;

}

void re_step(double * step, double * re_lsf, bool * mask, double * deltat, 
	double * xpr, double * xpl, double * ypf, double * ypb, double * zpu, double * zpd, 
	int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
{

	for(int pge_idx = 0; pge_idx < pges; ++pge_idx){
		for(int col_idx = 0; col_idx < cols; ++col_idx){
			for(int row_idx = 0; row_idx < rows; ++row_idx){
				
				int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);	

				int right = sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
				int right2 = sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
				int left = sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
				int left2 = sub2ind(row_idx, col_idx-2, pge_idx, rows, cols, pges);
				
				double xR, xL;
				eno_derivative(re_lsf[left2], re_lsf[left], re_lsf[ind], re_lsf[right], 
					re_lsf[right2], xpr[ind], xpl[ind], dx, &xR, &xL);

				int front = sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
				int front2 = sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
				int back = sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
				int back2 = sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);

				double yF, yB;
				eno_derivative(re_lsf[back2], re_lsf[back], re_lsf[ind], re_lsf[front],
					re_lsf[front2], ypf[ind], ypb[ind], dy, &yF, &yB);

				int up = sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
				int up2 = sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
				int down = sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
				int down2 = sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);

				double zU, zD;
				eno_derivative(re_lsf[down2], re_lsf[down], re_lsf[ind], re_lsf[up],
					re_lsf[up2], zpu[ind], zpd[ind], dz, &zU, &zD);

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