function RuncuSOLVERDemo
% Filename: RuncuSOLVERDemo.m
% Description: This function compiles the file cuSOLVERDemo.cu
% and compares the execution time of the solution of a system of
% linear equations Ax = b for x using MATLAB on the CPU and using
% cuSOLVER and cuBLAS on the GPU through the MEX file
% Authors: Ploskas, N., & Samaras, N.
% Syntax: RuncuSOLVERDemo
% Input: no inputs
% Output: the execution time of the computation of the solution of
% a system of linear equations Ax = b for x using MATLAB on the CPU
% and using cuSOLVER and cuBLAS on the GPU through the MEX file

% mexcuda cuSOLVERDemo.cu -lcublas -lcusolver
disp([	'Execution of the computation of the solution of a ' ...
	  	'system of linear equations Ax = b for x using MATLAB ' ...
		'on the CPU']);

A = rand(5000, 5000);
b = rand(5000, 1);

tic;
[~, n] = size(A);
[Q, R, P] = qr(A);
c = Q' * b;
y = R(1:n, 1:n) \ c(1:n);
x1 = P * y;
toc

disp([	'Execution of the computation of the solution of a ' ...
		'system of linear equations Ax = b for x using cuSOLVER ' ...
		'and cuBLAS on the GPU']);
gd = gpuDevice(); tic; x2 = cuSOLVERDemo(A, b); wait(gd); toc
