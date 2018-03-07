function RuncuSPARSEDemo
% Filename: RuncuSPARSEDemo.m
% Description: This function compiles the file cuSPARSEDemo.cu
% and compares the execution time of a sparse matrix-vector
% multiplication of a 5,000x5,000 array and a 5,000x1 vector
% using MATLAB on the CPU and using cuSPARSE on the GPU
% through the mex file
% Authors: Ploskas, N., & Samaras, N.
% Syntax: RuncuSPARSEDemo
% Input: no inputs
% Output: the execution time of a sparse matrix-vector
% multiplication of a 5,000x5,000 array and a 5,000x1 vector
% using MATLAB on the CPU and using cuSPARSE on the GPU
% through the mex file

% mexcuda cuSPARSEDemo.cu -lcusparse 

A = sprand(5000, 5000, 0.1);
b = rand(5000, 1);
disp([	'Execution of a sparse matrix-vector ' ...
		'multiplication of a 5,000x5,000 array and a ' ...
		'5,000x1 vector using MATLAB on the CPU']);

tic; x1 = A * b; toc

disp([	'Execution of a sparse matrix-vector ' ...
		'multiplication of a 5,000x5,000 array and a ' ...
		'5,000x1 vector using cuSPARSE on the GPU']);
gd = gpuDevice(); tic; x2 = cuSPARSEDemo(A, b); wait(gd); toc
