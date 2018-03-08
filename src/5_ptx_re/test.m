% generate -ptx code 
system('nvcc -ptx boundary_correction.cu');

% test mexcuda reinitialization scheme

addpath(genpath('..'))

% create a 3D level set function
xv = linspace(-5,5,64);
yv = xv;
zv = xv;

dx = xv(2) - xv(1);
dy = yv(2) - yv(1);
dz = zv(2) - zv(1);

[x, y, z] = meshgrid(xv,yv,zv);

fun = @(x,y,z) (0.1+(x-3.5).^2+(sqrt(y.^2+z.^2)-2).^2) .* (sqrt(x.^2/4+(z.^2+y.^2)/9)-1);

F = fun(x,y,z);

[rows, cols, pges] = size(F);
ThreadBlockSize = [rows, 4, 1];
num_ele = prod(size(F));

% define kernel function for boundary correction
k_bc = parallel.gpu.CUDAKernel('boundary_correction.ptx', 'boundary_correction.cu','boundary_correction');
k_bc.ThreadBlockSize = ThreadBlockSize;
k_bc.GridSize = [ceil(rows/ThreadBlockSize(1)), ...
				 ceil(cols/ThreadBlockSize(2)), ...
				 ceil(pges/ThreadBlockSize(3))];

% define kernel function for reinitialization step

re_step = parallel.gpu.CUDAKernel('boundary_correction.ptx', 'boundary_correction.cu','re_step');
re_step.ThreadBlockSize = ThreadBlockSize;
re_step.GridSize = [ceil(rows/ThreadBlockSize(1)), ...
				 	ceil(cols/ThreadBlockSize(2)), ...
					ceil(pges/ThreadBlockSize(3))];

function new_lsf = reini(k_bc, re_step, F, dx, dy, dz, rows, cols, pges, num_ele)

	Fgpu = gpuArray(F);

	mask = Fgpu<0;
	deltat = zeros(size(Fgpu),'gpuArray');

	xpr = ones(size(Fgpu), 'gpuArray') * dx;
	xpl = ones(size(Fgpu), 'gpuArray') * dx;
	ypf = ones(size(Fgpu), 'gpuArray') * dy;
	ypb = ones(size(Fgpu), 'gpuArray') * dy;
	zpu = ones(size(Fgpu), 'gpuArray') * dz;
	zpd = ones(size(Fgpu), 'gpuArray') * dz;

	[xpr, xpl, ypf, ypb, zpu, zpd] = feval(k_bc, xpr, xpl, ypf, ypb, zpu, zpd, Fgpu, num_ele, rows, cols, pges, dx, dy, dz); 

	minx = min(xpr,xpl);
	miny = min(ypf,ypb);
	minz = min(zpu,zpd);
	deltat = 0.3 * min(minx, min(miny,minz));

	step = zeros(size(Fgpu),'gpuArray');

	for i=1:100
		step = feval(re_step, step, Fgpu, mask, deltat, xpr, xpl, ypf, ypb, zpu, zpd, rows, cols, pges, dx, dy, dz, num_ele); 
		Ftmp = Fgpu - step;
		step = feval(re_step, step, Ftmp, mask, deltat, xpr, xpl, ypf, ypb, zpu, zpd, rows, cols, pges, dx, dy, dz, num_ele); 
		Fgpu = (Fgpu + Ftmp - step) / 2;
	end

	new_lsf = gather(Fgpu);

end




