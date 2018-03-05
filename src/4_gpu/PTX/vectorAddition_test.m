% example from : GPU Programming in MATLAB CH7.3
% test vectorAddtion.cu

% system('nvcc -ptx vectorAddition.cu');
k = parallel.gpu.CUDAKernel('vectorAddition.ptx', 'vectorAddition.cu');
n =10^7;
k.ThreadBlockSize = [k.MaxThreadsPerBlock, 1, 1];
k.GridSize = [ceil(n / k.ThreadBlockSize(1)), 1];

A = rand(n, 1);
B = rand(n, 1);

tic; C = A + B; toc

C = zeros(n, 1);

gd = gpuDevice(); tic; C = feval(k, A, B, C, n); wait(gd); toc

gpuA = gpuArray(A);
gpuB = gpuArray(B);
C = zeros(n, 1, 'gpuArray');

gd = gpuDevice(); tic; C = feval(k, A, B, C, n); wait(gd); toc
