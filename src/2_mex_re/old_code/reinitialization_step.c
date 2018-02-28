#include "mex.h"
#include "mexReinitialization.h"
#include "math.h"

void reinitialization_step(double *step, double *lsf, const boundary_modification* boundary, 
	const shift_mat * shift, const int num_ele, const int dx, const int dy, const int dz)
{
	//mexPrintf("now calculate the reinitialization step \n");

	double xR, xL, yF, yB, zU, zD;
	double f0, fr, fl, p2, p2r, p2l, p2m, sp;
	for (int i = 0, xjr = 0, xjl = 0, yjf = 0, yjb = 0, zju = 0, zjd = 0; i < num_ele; ++i)
	{	
		
		f0 = lsf[i];
		
		// compute xR & xL
		fr = lsf[shift->soXo[i]];
		fl = lsf[shift->soxo[i]];
		p2 = fr + fl - 2.0 * f0; // 2nd difference at current node	
		
		p2r = f0 + lsf[shift->soXo[shift->soXo[i]]] - 2.0 * fr; // 2nd difference at forward node
		p2m = 0.5 * MinMod(p2, p2r) / pow(dx,2); // 2nd derivative
		
		if (boundary->mxr[i]) {
			sp = boundary->xpr[xjr++];
			xR = (0-f0)/sp - sp * p2m; 
		}else{			
			xR = (fr-f0)/dx - dx * p2m;
		}

		p2l = f0 + lsf[shift->soxo[shift->soxo[i]]] - 2.0 * fl; // 2nd difference at backward node
		p2m = 0.5 * MinMod(p2, p2l) / pow(dx,2);
		if (boundary->mxl[i]) {
			sp = boundary->xpl[xjl++];
			xL = (f0-0)/sp + sp	* p2m;
		}else{
			xL = (f0-fl)/dx + dx * p2m;
		}
		// end computing xR & xL

		// compute yF & yB
		fr = lsf[shift->sYoo[i]];
		fl = lsf[shift->syoo[i]];
		p2 = fr + fl - 2.0 * f0; // 2nd difference at current node	
		
		p2r = f0 + lsf[shift->sYoo[shift->sYoo[i]]] - 2.0 * fr; // 2nd difference at forward node
		p2m = 0.5 * MinMod(p2, p2r) / pow(dy,2); // 2nd derivative
		
		if (boundary->myf[i]) {
			sp = boundary->ypf[yjf++];
			yF = (0-f0)/sp - sp * p2m; 
		}else{			
			yF = (fr-f0)/dy - dy * p2m;
		}

		p2l = f0 + lsf[shift->syoo[shift->syoo[i]]] - 2.0 * fl; // 2nd difference at backward node
		p2m = 0.5 * MinMod(p2, p2l) / pow(dy,2);
		if (boundary->myb[i]) {
			sp = boundary->ypb[yjb++];
			yB = (f0-0)/sp + sp	* p2m;
		}else{
			yB = (f0-fl)/dy + dy * p2m;
		}
		// end computing yF & yB

		// compute zU & zD
		fr = lsf[shift->sooZ[i]];
		fl = lsf[shift->sooz[i]];
		p2 = fr + fl - 2.0 * f0; // 2nd difference at current node	
		
		p2r = f0 + lsf[shift->sooZ[shift->sooZ[i]]] - 2.0 * fr; // 2nd difference at forward node
		p2m = 0.5 * MinMod(p2, p2r) / pow(dz,2); // 2nd derivative
		
		if (boundary->mzu[i]) {
			sp = boundary->zpu[zju++];
			zU = (0-f0)/sp - sp * p2m; 
		}else{			
			zU = (fr-f0)/dz - dz * p2m;
		}

		p2l = f0 + lsf[shift->sooz[shift->sooz[i]]] - 2.0 * fl; // 2nd difference at backward node
		p2m = 0.5 * MinMod(p2, p2l) / pow(dy,2);
		if (boundary->mzd[i]) {
			sp = boundary->zpd[zjd++];
			zD = (f0-0)/sp + sp	* p2m;
		}else{
			zD = (f0-fl)/dz + dz * p2m;
		}
		// end computing yF & yB

		// now calculate the time step
		// if inside
		if (boundary->mask[i]) {
			step[i] = ( sqrt(	Max2(pow(Min2(0,xL),2),pow(Max2(0,xR),2)) + 
								Max2(pow(Min2(0,yB),2),pow(Max2(0,yF),2)) + 
								Max2(pow(Min2(0,zD),2),pow(Max2(0,zU),2)) ) - 1)
						* boundary->deltat[i] * (-1.);
		} else{
			step[i] = ( sqrt(	Max2(pow(Max2(0,xL),2),pow(Min2(0,xR),2)) + 
								Max2(pow(Max2(0,yB),2),pow(Min2(0,yF),2)) + 
								Max2(pow(Max2(0,zD),2),pow(Min2(0,zU),2)) ) - 1)
						* boundary->deltat[i] * (1.);
		}

	}


}