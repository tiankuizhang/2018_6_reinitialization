% reinitialization scheme with PTX gpu implementation
% for a 64x64x64 array, 0.04s is used: faster than the mexcuda scheme!
function new_lsf = ptx_re(k_bc, re_step, F, dx, dy, dz, rows, cols, pges, num_ele)

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




