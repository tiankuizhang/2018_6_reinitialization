// computing routine for reinitialization
#include "mex.h"
#include "mexReinitialization.h"
#include "math.h"


#define Discriminant(p2,v0,v2) ( pow((0.5*p2-v0-v2),2) - 4.*v0*v2 )
#define dist(disc,ds,p2,v0,v2) ( ds * (0.5 + (v0-v2 - sign(v0-v2)*sqrt(disc)) / p2 ) )
#define dist_turn(ds,v0,v2) ( ds * v0 / (v0 - v2) )
#define epsilon 10e-10

void Reinitialization(double *re_lsf, const double *lsf, const shift_mat * shift, 
	const int rows, const int cols, const int pages, const int num_ele, const int dx, const int dy, const int dz)
{

	// binary mask. ture if there is boundary to the right
	bool *mask = (bool *)mxCalloc((mwSize)num_ele,sizeof(bool)); // ture for inside, false outside
	double *deltat = (double *)mxCalloc((mwSize)num_ele,sizeof(double)); // time step
	bool *mxr = (bool *)mxCalloc((mwSize)num_ele,sizeof(bool));
	bool *mxl = (bool *)mxCalloc((mwSize)num_ele,sizeof(bool));
	bool *myf = (bool *)mxCalloc((mwSize)num_ele,sizeof(bool));
	bool *myb = (bool *)mxCalloc((mwSize)num_ele,sizeof(bool));
	bool *mzu = (bool *)mxCalloc((mwSize)num_ele,sizeof(bool));
	bool *mzd = (bool *)mxCalloc((mwSize)num_ele,sizeof(bool));
	int num_xp = 0; // number of pair of nodes crossing boundary
	int num_yf = 0;
	int num_zu = 0;
	double min_ds = Min2(Min2(dx,dy),dz);
	for (int i = 0; i < num_ele; ++i){
		re_lsf[i] = lsf[i];
		mask[i] = (lsf[i]<0) ;
		deltat[i] = min_ds;
		mxr[i] = (lsf[i]*lsf[shift->soXo[i]] < 0); 
		mxl[i] = (lsf[i]*lsf[shift->soxo[i]] < 0);
		myf[i] = (lsf[i]*lsf[shift->sYoo[i]] < 0);
		myb[i] = (lsf[i]*lsf[shift->syoo[i]] < 0);
		mzu[i] = (lsf[i]*lsf[shift->sooZ[i]] < 0);
		mzd[i] = (lsf[i]*lsf[shift->sooz[i]] < 0);
		num_xp += (int)(mxr[i]);
		num_yf += (int)(myf[i]);
		num_zu += (int)(mzu[i]);
	}
	
	double *xpr = (double *)mxCalloc((mwSize)num_xp,sizeof(double)); // distance from boundary : xpr
	double *xpl = (double *)mxCalloc((mwSize)num_xp,sizeof(double));
	double *ypf = (double *)mxCalloc((mwSize)num_yf,sizeof(double));
	double *ypb = (double *)mxCalloc((mwSize)num_yf,sizeof(double));
	double *zpu = (double *)mxCalloc((mwSize)num_zu,sizeof(double));
	double *zpd = (double *)mxCalloc((mwSize)num_zu,sizeof(double));
	for (int i = 0, xj = 0, yj = 0, zj = 0; i < num_ele; ++i){
		double f0, f2, p2;
		// boundary modification in x direction
		if (mxr[i]) {
			f0 = lsf[i]; // grab the left node value
			f2 = lsf[shift->soXo[i]]; // grab the right node value
			double p2xl = lsf[shift->soxo[i]] + f2 - 2.0 * f0; // 2nd difference on the left node
			double p2xr = f0 + lsf[shift->soXo[shift->soXo[i]]] - 2.0 * f2; // 2nd difference right
			p2 = MinMod(p2xl, p2xr);
			if (p2>epsilon){
				//double disc = Discriminant(p2x, f0, f2);
				xpr[xj] = dist(Discriminant(p2,f0,f2),dx,p2,f0,f2);
			} else{
				xpr[xj] = dist_turn(dx,f0,f2);
			}
			xpl[xj] = dx - xpr[xj];
			deltat[i] = Min2(deltat[i],xpr[xj]); // choose the smallest time step
			deltat[shift->soXo[i]] = Min2(deltat[shift->soXo[i]],xpl[xj]);
			++xj;
		} 
		// boundary modification in y dirction
		if (myf[i]) {
			f0 = lsf[i]; // grab the back node value
			f2 = lsf[shift->sYoo[i]]; // grab the front node value
			double p2yb = lsf[shift->syoo[i]] + f2 - 2.0 * f0; // 2nd difference on the back node
			double p2yf = f0 + lsf[shift->sYoo[shift->sYoo[i]]] - 2.0 * f2; // 2nd difference on the front node
			p2 = MinMod(p2yb, p2yf);
			if (p2>epsilon){
				ypf[yj] = dist(Discriminant(p2,f0,f2),dy,p2,f0,f2);
			} else{
				ypf[yj] = dist_turn(dy,f0,f2);
			}
			ypb[yj] = dy - ypf[yj];
			deltat[i] = Min2(deltat[i],ypf[yj]);
			deltat[shift->sYoo[i]] = Min2(deltat[shift->sYoo[i]],ypb[yj]);
			++yj;
		}
		// boundary modification in z direction
		if (mzu[i]) {
			f0 = lsf[i]; // grab the down node value
			f2 = lsf[shift->sooZ[i]]; // grab the up node value
			double p2zd = lsf[shift->sooz[i]] + f2 - 2.0 * f0; // 2nd difference on the down node
			double p2zu = f0 + lsf[shift->sooZ[shift->sooZ[i]]] - 2.0 * f2; // 2nd diference on the up node
			p2 = MinMod(p2zd, p2zu);
			if (p2>epsilon){
				zpu[zj] = dist(Discriminant(p2,f0,f2),dz,p2,f0,f2);
			} else{
				zpu[zj] = dist_turn(dz,f0,f2);
			}
			zpd[zj] = dz - zpu[zj];
			deltat[i] = Min2(deltat[i],zpu[zj]);
			deltat[shift->sooZ[i]] = Min2(deltat[shift->sooZ[i]],zpd[zj]);
			++zj;
		}
	}

	for (int i = 0; i < num_ele; ++i)
	{
		deltat[i] = 0.3 * deltat[i];
	}

	// now gather all loop invariants
	boundary_modification boundary;
	boundary.mask = mask;
	boundary.mxr = mxr;
	boundary.mxl = mxl;
	boundary.myf = myf;
	boundary.myb = myb;
	boundary.mzu = mzu;
	boundary.mzd = mzd;
	boundary.xpr = xpr;
	boundary.xpl = xpl;
	boundary.ypf = ypf;
	boundary.ypb = ypb;
	boundary.zpu = zpu;
	boundary.zpd = zpd;
	boundary.deltat = deltat;

	// iteraction with runge kutta to update the distance map
	// re_lsf will serve as the current distance map
	double *step = (double *)mxCalloc((mwSize)num_ele, sizeof(double)); // hold time step
	double *intermediate = (double *)mxCalloc((mwSize)num_ele, sizeof(double)); // hold intermediate step

	//reinitialization_step(step, re_lsf, &boundary, shift, num_ele, dx, dy, dz);
	for (int i = 0; i < 100; ++i)
	{
		/* code */
		reinitialization_step(step, re_lsf, &boundary, shift, num_ele, dx, dy, dz);
		for(int j = 0; j < num_ele; ++j)
			intermediate[j] = re_lsf[j] - step[j];
		reinitialization_step(step, intermediate, &boundary, shift, num_ele, dx, dy, dz);
		for(int j = 0; j < num_ele; ++j)
			re_lsf[j] = (re_lsf[j] + intermediate[j] - step[j]) / 2;
	}



	// end of iteration

	mxFree(mask);	
	mxFree(deltat);
	mxFree(mxr);
	mxFree(mxl);
	mxFree(myf);
	mxFree(myb);
	mxFree(mzu);
	mxFree(mzd);
	mxFree(xpr);
	mxFree(xpl);
	mxFree(ypf);
	mxFree(ypb);
	mxFree(zpu);
	mxFree(zpd);
	mxFree(step);
	mxFree(intermediate);

	//mexPrintf("%u pairs of nodes in x dirction cross boundary \n", num_xp);
	//mexPrintf("%u pairs of nodes in y dirction cross boundary \n", num_yf);
	//mexPrintf("%u pairs of nodes in z dirction cross boundary \n", num_zu);
}