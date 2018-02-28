// mexReinitialization.c
// implement the reinitialization scheme in 3D
// in hope of reducing computation time 

#include "mex.h"
#include "mexReinitialization.h"

// enumerate input index
enum Input {
	level_set_function = 0,
	shift_matrix = 1,
	grid_spacing = 2
};

// enumerate output index
enum out_put{
	reinitialized_lsf = 0
};

// enumerate the shift matrix struct field name/index
enum field_name{
	soXo = 0,
	soxo = 1,
	sYoo = 2,
	syoo = 3,
	sooZ = 4,
	sooz = 5
};

// the gateway function
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{	
	mxClassID category;
	
	if(nrhs != 3){
		mexErrMsgIdAndTxt("mexReinitialization:wrong_number_of_inputs",
			"expecting 3 inputs");
	}
	
	// assign level set function to LSF
	double *lsf;

	mwSize number_of_dimensions;
	const mwSize *dimension_array;
	size_t number_of_elements_lsf;

	category = mxGetClassID(prhs[level_set_function]);
	number_of_dimensions = mxGetNumberOfDimensions(prhs[level_set_function]);
	dimension_array = mxGetDimensions(prhs[level_set_function]);
	number_of_elements_lsf = mxGetNumberOfElements(prhs[level_set_function]);
	if (category != mxDOUBLE_CLASS || number_of_dimensions != (mwSize)3){
		mexErrMsgIdAndTxt("mexReinitialization:Invalid_Input",
			"Argument %d must be a 3 dimension array of double precision!",
			level_set_function);
	}
	lsf = (double *)mxGetData(prhs[level_set_function]);
	// finish assigning LSF
	
	// assign shift_mat to the struct shift in C
	shift_mat shift;
	
	mwSize total_number_of_elements;
	int number_of_fields;

	category = mxGetClassID(prhs[shift_matrix]);
	total_number_of_elements = mxGetNumberOfElements(prhs[shift_matrix]);
	number_of_fields = mxGetNumberOfFields(prhs[shift_matrix]);

	if (category != mxSTRUCT_CLASS || total_number_of_elements != (mwSize)1 || number_of_fields != (int)6){
		mexErrMsgIdAndTxt("mexReinitialization:Invalid_Input",
			"Argument %d must be a struct containing 6 index matrices shifted in X,x,Y,y,Z,z directions!", 
			shift_matrix);
	}
	shift.soXo = (int *)mxGetData(mxGetFieldByNumber(prhs[shift_matrix], (mwIndex)0, soXo));
	shift.soxo = (int *)mxGetData(mxGetFieldByNumber(prhs[shift_matrix], (mwIndex)0, soxo));
	shift.sYoo = (int *)mxGetData(mxGetFieldByNumber(prhs[shift_matrix], (mwIndex)0, sYoo));
	shift.syoo = (int *)mxGetData(mxGetFieldByNumber(prhs[shift_matrix], (mwIndex)0, syoo));
	shift.sooZ = (int *)mxGetData(mxGetFieldByNumber(prhs[shift_matrix], (mwIndex)0, sooZ));
	shift.sooz = (int *)mxGetData(mxGetFieldByNumber(prhs[shift_matrix], (mwIndex)0, sooz));	
	// finish assigning shift_mat

	// assign grid spacing array
	double *ds;
	size_t rows, cols;

	category = mxGetClassID(prhs[grid_spacing]);
	rows = mxGetM(prhs[grid_spacing]);
	cols = mxGetN(prhs[grid_spacing]);
	if (category != mxDOUBLE_CLASS || rows != (size_t)1 || cols != (size_t)3){
		mexErrMsgIdAndTxt("mexReinitialization:Invalid_Input",
			"Argument %d must be a 1X3 double array of the grid spacing",
			grid_spacing);
	}
	ds = (double *)mxGetData(prhs[grid_spacing]);
	// finish assigning spacing array
	

	// create an output array
	double *re_lsf; // pointer to reinitialized level set function
	plhs[reinitialized_lsf] = mxCreateNumericArray(number_of_dimensions, dimension_array, mxDOUBLE_CLASS, mxREAL);
	re_lsf = (double *)mxGetData(plhs[reinitialized_lsf]);
	//  finish creating output array

	//Reinitialization();
	//mexPrintf("number of elements :%d \n", number_of_elements_lsf);
	//mexPrintf("dx:%f dy:%f dz:%f \n",ds[0],ds[1],ds[2]);

	// computation routine
	Reinitialization(re_lsf, lsf, &shift, dimension_array[0], dimension_array[1], dimension_array[2], 
		number_of_elements_lsf, ds[0], ds[1], ds[2]);

}