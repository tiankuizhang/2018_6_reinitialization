#ifndef _mexReinitialization_H
#define _mexReinitialization_H

#define MinMod(x,y) ( (x*y<0) ? 0 : (fabs(x)<fabs(y) ? x : y ) ) 
#define Min2(x,y) ( (x<y) ? x : y )
#define Max2(x,y) ( (x<y) ? y : x )
#define sign(x) ( x>0 ? 1. : -1. )

// struct defined to avoid passing a long list of pointers to the reinitialization scheme
typedef struct
{	
	int *soXo;
	int *soxo;
	int *sYoo;
	int *syoo;
	int *sooZ;
	int *sooz;
} shift_mat;

// struct defined to avoid passing a long list of pointers to the step function
typedef struct {
	bool *mask;
	bool *mxr;
	bool *mxl;
	bool *myf;
	bool *myb;
	bool *mzu;
	bool *mzd;
	double *xpr;
	double *xpl;
	double *ypf;
	double *ypb;
	double *zpu;
	double *zpd;
	double *deltat;
} boundary_modification; 

// prototype for the Reinitialization function doing the actual calculation
// re_lsf		reinitialized level set function
// lsf 			initial level set function 
// shift 		index shift matrix
// r,c,p 		matrix dimension array
// num_ele 		total number of elements
// dx,dy,dz 	grid spacing
void Reinitialization(double *re_lsf, const double *lsf, const shift_mat * shift, 
	const int rows, const int cols, const int pages, const int num_ele, const int dx, const int dy, const int dz);


// calculate one reinitialization step 
void reinitialization_step(double *  step, double *  lsf, const boundary_modification*  boundary,
	const shift_mat *  shift , const int num_ele, const int dx, const int dy, const int dz);

#endif


