function RuncuRANDDemo
% Filename: RuncuRANDDemo.m
% Description: This function compiles the file cuRANDDemo.cu
% and compares the execution time of the computation of the area
% of the definite integral int(x ∗ (x - 2) ^ 6, x = 0..2) using
% MATLAB on the CPU and using cuRAND on the GPU through the MEX
% file
% Authors: Ploskas, N., & Samaras, N.
% Syntax: RuncuRANDDemo
% Input: no inputs
% Output: the execution time of the computation of the area of
% the definite integral int(x ∗ (x - 2) ^ 6, x = 0..2) using
% MATLAB on the CPU and using cuRAND on the GPU through the MEX
% file

% mexcuda cuRANDDemo.cu -lcurand 

disp([	'Execution of the computation of the area of the ' ...
		'the definite integral int(x ∗ (x - 2) ^ 6, x = 0..2) ' ...
		'using MATLAB on the CPU']);

tic;
n = 1000000;
x = 2 * rand(n, 1);
y = 8 * rand(n, 1);
counter = sum(y < (x .* (x - 2) .^ 6));
area = counter / n * 2 * 8;
toc

disp([	'Execution of the computation of the area of the ' ...
		'the definite integral int(x ∗ (x - 2) ^ 6, x = 0..2) ' ...
		'using cuRAND on the GPU']);

gd = gpuDevice(); tic; area = cuRANDDemo(); wait(gd); toc

