
function RunNPPDemo(filename)
% Filename: RunNPPDemo.m
% Description: This function compiles the file NPPDemo.cu
% and compares the execution time of a 2-D box filteringo
% on a 3,867x5,789 pgm image using MATLAB on the CPU and
% using NPP on the GPU through the MEX file
% Authors: Ploskas, N., & Samaras, N.
% Syntax: RunNPPDemo
% Input:
% -- filename: the full path a pgm image
% Output: the execution time of a 2-D box filtering on a
% 3,867x5,789 pgm image using MATLAB on the CPU and
% using NPP on the GPU through the MEX file

mexcuda NPPDemo.cu -lnppi -lfreeimage ...
	-I/cm/shared/apps/cuda80/sdk/8.0.61/7_CUDALibraries/common/UtilNPP ...
    -I/cm/shared/apps/cuda80/sdk/8.0.61/7_CUDALibraries/common/FreeImage/include ...
    -L/cm/shared/apps/cuda80/sdk/8.0.61/7_CUDALibraries/common/FreeImage/lib/linux/x86_64
			

disp([	'Execution of a 2-D box filtering on a ' ...
		'3,867x5,789 pgm image using MATLAB on the CPU']);

tic;
A = imread(filename);
B = imboxfilt(A, 11, 'padding', 'replicate');
[pathstr, name, ext] = fileparts(filename);
newFilename = strcat(pathstr, '\', name, '_boxFilterMATLAB', ext);
imwrite(B, newFilename);
toc

disp([	'Execution of a 2-D box filtering on a ' ...
		'3,867x5,789 pgm image using NPP on the GPU through']);

gd = gpuDevice(); tic; NPPDemo(filename); wait(gd); toc
