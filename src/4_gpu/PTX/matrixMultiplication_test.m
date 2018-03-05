% example from GPU computing with MATLAB
% system('nvcc -ptx matrixMultiplication.cu');

k = parallel.gpu.CUDAKernel('matrixMultiplication.ptx', 'matrixMultiplication.cu');
n = 5000;
k.ThreadBlockSize = [10, 10, 1];
k.GridSize = [100, 100];
A = rand(n);
B = rand(n);

tic; C = A * B; toc

C = zeros(n);

gd = gpuDevice(); tic; C = feval(k, A, B, C, n); wait(gd); toc

gpuA = gpuArray(A);
gpuB = gpuArray(B);
C = zeros(n, 'gpuArray');
gpd = gpuDevice(); tic; C = feval(k, gpuA, gpuB, C, n); wait(gd); toc
