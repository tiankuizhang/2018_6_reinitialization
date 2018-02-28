#ifndef _MEXRE_H
#define _MEXRE_H

inline int sub2ind(int row_idx, int col_idx, int pge_idx, int rows, int cols, int pges);
inline double max2(double x, double y);
inline double min2(double x, double y);
inline double min_mod(double x, double y);
inline double sign(double x);
inline double sp(double v1, double v0, double v2, double v3, double ds);
inline double Discriminant(double p2, double v0, double v2);
inline double dist(double disc, double ds,double p2,double v0, double v2);
inline double dist_turn(double ds, double v0, double v2);

void reinitialization(double * re_lsf, double * lsf, int rows, int cols, int pges, 
	double dx, double dy, double dz, int num_ele);

void re_step(double * step, double * re_lsf, bool * mask, double * deltat, 
	double * xpr, double * xpl, double * ypf, double * ypb, double * zpu, double * zpd, 
	int rows, int cols, int pges, double dx, double dy, double dz, int num_ele);

#endif