function RuncuFFTDemo
% Filename: RuncuFFTDemo.m
% Description: This function compiles the file cuFFTDemo.cu
% and compares the execution time of a two-dimensional
% discrete Fourier transform of a 5,000x5,000 array using
% MATLAB on the CPU and using cuFFT on the GPU through 
% the MEX file
% Authors: Ploskas, N., & Samaras, N.
% Syntax: RuncuFFTDemo
% Input: no inputs
% Output: the execution time of a two-dimensional discrete
% Fourier transform of a 5,000x5,000 array using MATLAB on 
% the CPU and using cuBLAS on the GPU through the MEX file

mexcuda cuFFTDemo.cu -lcufft

A = rand(5000);
disp(['Execution of a two-dimensional discrete' ...
	  'Fourier transform of a 5,000x5,000 array' ...
	  'using MATLAB on the CPU']);

tic; B1 = fft2(A); toc
disp(['Execution of a two-dimensional discrete' ...
	  'Fourier transform of a 5,000x5,000 array' ...
	  'using cuFFT on the GPU']);
gd = gpuDevice(); tic; B2 = cuFFTDemo(A);wait(gd);toc
