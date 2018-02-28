#include "mex.h"
#include "mexRe.h"
#include "math.h"


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