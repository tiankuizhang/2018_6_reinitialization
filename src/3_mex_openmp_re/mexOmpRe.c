/* implementation of reinitialization scheme
 * trying to use openmp 
 * calling syntax : new_lsf = mexRe(lsf,dx,dy,dz)
 * lsf : old 3d level set function
 * dx,dy,dz : grid spacing in x,y,z direction
 * 
 * (failed to set CFLAGS)compile on linux : mex mexOmpRe.c Re.c CXXFLAGS="$CXXFLAGS -std=c99"
 * (worked)compile on linux : mex -v mexOmpRe.c Re.c re_step.c func.c CFLAGS="$CFLAGS -std=c99 -O3"
 */

#include "mex.h"
#include "omp.h"
#include "mexOmpRe.h"

enum in_put {
	level_set_function = 0,
	grid_spacing_x = 1,
	grid_spacing_y = 2,
	grid_spacing_z = 3
};

enum out_put {
	reinitialized_lsf = 0
};

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{	

	mxClassID category;

	if(nrhs != 4){
		mexErrMsgIdAndTxt("mexReinitialization:wrong_number_of_inputs",
			"expecting 4 inputs");
	}

	// assign level set function to LSF
	double *lsf;

	mwSize number_of_dimensions;
	const mwSize *dimension_array;
	size_t num_ele;

	category = mxGetClassID(prhs[level_set_function]);
	number_of_dimensions = mxGetNumberOfDimensions(prhs[level_set_function]);
	dimension_array = mxGetDimensions(prhs[level_set_function]);
	num_ele = mxGetNumberOfElements(prhs[level_set_function]);
	if (category != mxDOUBLE_CLASS || number_of_dimensions != (mwSize)3){
		mexErrMsgIdAndTxt("mexReinitialization:Invalid_Input",
			"Argument %d must be a 3 dimension array of double precision!",
			level_set_function);
	}
	lsf = (double *)mxGetData(prhs[level_set_function]);
	// finish assigning LSF

	// assign size of level set function
	int cols = dimension_array[0];
	int rows = dimension_array[1];
	int pges = dimension_array[2];

	mexPrintf("size of level_set_function: (%d,%d,%d)\n", cols, rows, pges);

	// assign grid spacing array
	double dx = mxGetScalar(prhs[grid_spacing_x]);
	double dy = mxGetScalar(prhs[grid_spacing_y]);
	double dz = mxGetScalar(prhs[grid_spacing_z]);

	mexPrintf("grid_spacing: (%f,%f,%f)\n", dx, dy, dz);

	// create an output array
	double *re_lsf; // pointer to reinitialized level set function
	plhs[reinitialized_lsf] = mxCreateNumericArray(number_of_dimensions, 
		dimension_array, mxDOUBLE_CLASS, mxREAL);
	re_lsf = (double *)mxGetData(plhs[reinitialized_lsf]);
	//  finish creating output array

	reinitialization(re_lsf, lsf, rows, cols, pges, dx, dy, dz, num_ele);

}
