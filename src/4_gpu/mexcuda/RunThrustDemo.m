function RunThrustDemo
% Filename: RunThrustDemo.m
% Description: This function compiles the file ThrustDemo.cu
% and compares the execution time of the calculation of
% the sum of a 100,000,000x1 vector elements using MATLAB
% on the CPU and using Thrust on the GPU through the MEX
% file
% Authors: Ploskas, N., & Samaras, N.
% Syntax: RunThrustDemo
% Input: no inputs
% Output: the execution time of the calculation of the
% sum of a 100,000,000x1 vector elements using MATLAB
% on the CPU and using Thrust on the GPU through the MEX
% file

mexcuda ThrustDemo.cu
A = rand(100000000, 1);
disp([	'Execution of the calculation of the sum of ' ...
		'a 100,000,000x1 vector elements using MATLAB ' ...
		'on the CPU']);

tic; b1 = sum(A); toc

disp([	'Execution of the calculation of the sum of ' ...
		'a 100,000,000x1 vector elements using Thrust ' ...
		'on the GPU']);

gd = gpuDevice(); 
tic; 
b2 = ThrustDemo(A); 
wait(gd); 
toc
