/*==========================================================
 * Filename: NPPDemo.cu
 * Description: This function implements a 2-D box filtering
 * of a pgm image using NPP library (MEX file that contains
 * CUDA code and takes as inputs MATLAB arrays)
 * Authors: Ploskas, N., & Samaras, N.
 * Syntax: NPPDemo(filename)
 * Input:
 * -- filename: the path to a pgm image
 * Output: a new pgm image is stored at the same location
 * with the original image with a name name_boxFilter.pgm,
 * where name is the name of the original image
 *========================================================*/

#include "mex.h"
#include <cuda_runtime.h>
#include <npp.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <ImagesCPU.h>
#include <ImagesNPP.h>
#include <ImageIO.h>
#include <Exceptions.h>
#include <helper_string.h>
#include <helper_cuda.h>

/*
 * The gateway function
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	
	/* declare all variables */
	char *filename;
	std::string name, sResultFilename;

	/* define error messages */
	char const * const errId = "parallel:gpu:NPPDemo:InvalidInput";
	char const * const errMsg1 = "Invalid input to MEX file.";
	char const * const errMsg2 = "An exception occurred.";

	/* check input data */
	if(nrhs != 1){
		mexErrMsgIdAndTxt(errId, errMsg1);
	}

	/* copy the image filename from prhs[0] into filename */
	filename = mxArrayToString(prhs[0]);
	name = filename;

	try{
		/* try to open the image */
		std::ifstream infile(name.data(), std::ifstream::in);
		if (!infile.good()){
			infile.close();
			mexErrMsgIdAndTxt(errId, errMsg1);
		}

		/* create the name of the final image */
		sResultFilename = name;
		std::string::size_type dot = sResultFilename.rfind('.');
		if (dot != std::string::npos){
			sResultFilename = sResultFilename.substr(0, dot);
		}
		sResultFilename += "_boxFilterNPP.pgm";

		/* declare a host image object for an 8-bit grayscale
		   image */
		npp::ImageCPU_8u_C1 oHostSrc;

		/* load the image from disk */
		npp::loadImage(name, oHostSrc);

		/* declare a device image and copy construct from the
		   host image */
		npp::ImageNPP_8u_C1 oDeviceSrc(oHostSrc);

		/* create struct with box-filter mask size */
		NppiSize oMaskSize = {11, 11};
		NppiSize oSrcSize = {(int)oDeviceSrc.width(), (int)oDeviceSrc.height()};
		NppiPoint oSrcOffset = {0, 0};

		/* create struct with ROI size */
		NppiSize oSizeROI = {(int)oDeviceSrc.width(), (int)oDeviceSrc.height()};
		
		/* allocate device image of appropriately reduced size */
		npp::ImageNPP_8u_C1 oDeviceDst(oSizeROI.width, oSizeROI.height);

		/* set anchor point inside the mask to (oMaskSize.width / 2,
		   oMaskSize.height / 2). It should round down when odd */
		NppiPoint oAnchor = {oMaskSize.width / 2, oMaskSize.height / 2};

		/* perform the box filter */
		NPP_CHECK_NPP(nppiFilterBoxBorder_8u_C1R(oDeviceSrc.data(),
				   	oDeviceSrc.pitch(),
				   	oSrcSize,
					oSrcOffset,
					oDeviceDst.data(),
					oDeviceDst.pitch(),
					oSizeROI,
					oMaskSize,
					oAnchor,
					NPP_BORDER_REPLICATE) );

		/* declare a host image for the result and copy the
		   device result data into it*/
		npp::ImageCPU_8u_C1 oHostDst(oDeviceDst.size());
		oDeviceDst.copyTo(oHostDst.data(), oHostDst.pitch());
		npp::saveImage(sResultFilename, oHostDst);

		/* destroy the data on the device */
		nppiFree(oDeviceSrc.data());
		nppiFree(oDeviceDst.data());
	}
	catch (npp::Exception &rException){
	   	mexErrMsgIdAndTxt(errId, errMsg2);
	}
	catch (...){
	   	mexErrMsgIdAndTxt(errId, errMsg2);
	}
}
