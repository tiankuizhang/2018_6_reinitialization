function reinitialization(obj, Distance)

	Epsilon = 10^(-10);
	Ny = obj.GD3.mrows;
	Nx = obj.GD3.ncols;
	Nz = obj.GD3.lshts;

	Dx = obj.GD3.Dx;
	Dy = obj.GD3.Dy;
	Dz = obj.GD3.Dz;

	xpr = Dx*ones(Ny,Nx,Nz); xpl = xpr; % grid spacing in x direction
	ypf = Dy*ones(Ny,Nx,Nz); ypb = ypf; % grid spacing in y direction
	zpu = Dz*ones(Ny,Nx,Nz); zpd = zpu; % grid spacing in z direction

	d0 = Distance;

%	dl0 = circshift(d0, [0 1 0]);  dl0(:,1,:)   = d0(:,1,:);	% x:left
%	dr0 = circshift(d0, [0 -1 0]); dr0(:,end,:) = d0(:,end,:);	% x:right
%	df0 = circshift(d0, [-1 0 0]); df0(end,:,:) = d0(end,:,:);	% y:front
%	db0 = circshift(d0, [1 0 0]);  db0(1,:,:)   = d0(1,:,:);	% y:back
%	du0 = circshift(d0, [0 0 -1]); du0(:,:,end) = d0(:,:,end);	% z:up
%	dd0 = circshift(d0, [0 0 1]);  dd0(:,:,1)   = d0(:,:,1);	% z:down

	dl0 = d0(obj.GD3.oxo); dl0(:,1,:)   = d0(:,1,:);	% x:left
	dr0 = d0(obj.GD3.oXo); dr0(:,end,:) = d0(:,end,:);	% x:right
	df0 = d0(obj.GD3.Yoo); df0(end,:,:) = d0(end,:,:);	% y:front
	db0 = d0(obj.GD3.yoo); db0(1,:,:)   = d0(1,:,:);	% y:back
	du0 = d0(obj.GD3.ooZ); du0(:,:,end) = d0(:,:,end);	% z:up
	dd0 = d0(obj.GD3.ooz); dd0(:,:,1)   = d0(:,:,1);	% z:down

		% x direction

		mxr = d0.*dr0<0; % boundary points to the right of the boundary in x direction
		mxl = circshift(mxr, [0 1 0]); % boundary points to the left of the boundary in x direction
			
		P2x0r = (dl0(mxr)+dr0(mxr)-2*d0(mxr));
		P2x0l = (dl0(mxl)+dr0(mxl)-2*d0(mxl));

		mmodPx = zeros(Ny,Nx,Nz); 
		mmodPx(mxr) = MinMod(P2x0r, P2x0l);			
		mmodPx(mxl) = mmodPx(mxr); %mmodPx(mxl) = MinMod(P2x0r, P2x0l);

		% boundary modification

		mxr1 = mxr & abs(mmodPx)>Epsilon; % not a turning point
		mxr2 = mxr & ~mxr1; % turning point

		t1 = mmodPx(mxr1); t2 = d0(mxr1); t3 = dr0(mxr1); % some scaffod for compuatation
		Disc = (0.5*t1-t2-t3).^2 - 4*t2.*t3; % Discriminat
		xpr(mxr1) = Dx * (0.5 + (t2-t3-sign(t2-t3).*sqrt(Disc)) ./t1 );
		xpr(mxr2) = Dx * d0(mxr2) ./ (d0(mxr2)-dr0(mxr2)); % distacne from grid mxr to the boundary point

		% boundary modification
		mxl1 = mxl & abs(mmodPx)>Epsilon; % not a turning point
		mxl2 = mxl & ~mxl1; % turning point

		t1 = mmodPx(mxl1); t2 = d0(mxl1); t3 = dl0(mxl1); % some scaffod for compuatation
		Disc = (0.5*t1-t2-t3).^2 - 4*t2.*t3; % Discriminat
		xpl(mxl1) = Dx * (0.5 + (t2-t3-sign(t2-t3).*sqrt(Disc)) ./t1 );
		xpl(mxl2) = Dx * d0(mxl2) ./ (d0(mxl2)-dl0(mxl2)); % distacne from grid mxr to the boundary point

	% y direction

		myf = d0.*df0<0; % boundary points to the back of the boundary in y direction
		myb = circshift(myf, [1 0 0]);
					
		P2y0f = (df0(myf)+db0(myf)-2*d0(myf));
		P2y0b = (df0(myb)+db0(myb)-2*d0(myb));

		mmodPy = zeros(Ny,Nx,Nz);  
		mmodPy(myf) = MinMod(P2y0f, P2y0b);
		mmodPy(myb) = mmodPy(myf); %mmodPy(myb) = MinMod(P2y0f, P2y0b);

		% boundary modification
		myf1 = myf & abs(mmodPy)>Epsilon; % not a turning point
		myf2 = myf & ~myf1; % turning point

		t1 = mmodPy(myf1); t2 = d0(myf1); t3 = df0(myf1); % some scaffod for compuatation
		Disc = (0.5*t1-t2-t3).^2 - 4*t2.*t3; % Discriminat
		ypf(myf1) = Dy * (0.5 + (t2-t3-sign(t2-t3).*sqrt(Disc)) ./t1 );
		ypf(myf2) = Dy * d0(myf2) ./ (d0(myf2)-df0(myf2)); % distacne from grid mxr to the boundary point

		% boundary modification
		myb1 = myb & abs(mmodPy)>Epsilon; % not a turning point
		myb2 = myb & ~myb1; % turning point

		t1 = mmodPy(myb1); t2 = d0(myb1); t3 = db0(myb1); % some scaffod for compuatation
		Disc = (0.5*t1-t2-t3).^2 - 4*t2.*t3; % Discriminat
		ypb(myb1) = Dy * (0.5 + (t2-t3-sign(t2-t3).*sqrt(Disc)) ./t1 );
		ypb(myb2) = Dy * d0(myb2) ./ (d0(myb2)-db0(myb2)); % distacne from grid mxr to the boundary point

	% z direction

		mzu = d0.*du0<0; % boundary points to the downward of the boundary in the z direction
		mzd = circshift(mzu, [0 0 1]);

		

		P2z0u = (du0(mzu)+dd0(mzu)-2*d0(mzu));
		P2z0d = (du0(mzd)+dd0(mzd)-2*d0(mzd));

		mmodPz = zeros(Ny,Nx,Nz);
		mmodPz(mzu) = MinMod(P2z0u,P2z0d);
		mmodPz(mzd) = mmodPz(mzu);

		% boundary modification
		mzu1 = mzu & abs(mmodPz)>Epsilon; % not a turning point
		mzu2 = mzu & ~mzu1; % turning point

		t1 = mmodPz(mzu1); t2 = d0(mzu1); t3 = du0(mzu1); % some scaffod for computation
		Disc = (0.5*t1-t2-t3).^2 - 4*t2.*t3; % Discriminat
		zpu(mzu1) = Dz * (0.5 + (t2-t3-sign(t2-t3).*sqrt(Disc))./t1);
		zpu(mzu2) = Dz * d0(mzu2) ./ (d0(mzu2)-du0(mzu2)); % distacne from grid mxr to the boundary point in z direction

		% boundary modification
		mzd1 = mzd & abs(mmodPz)>Epsilon; % not a turning point
		mzd2 = mzd & ~mzd1; % turning point

		t1 = mmodPz(mzd1); t2 = d0(mzd1); t3 = dd0(mzd1); % some scaffod for computation
		Disc = (0.5*t1-t2-t3).^2 - 4*t2.*t3; % Discriminat
		zpd(mzd1) = Dz * (0.5 + (t2-t3-sign(t2-t3).*sqrt(Disc))./t1);
		zpd(mzd2) = Dz * d0(mzd2) ./ (d0(mzd2)-dd0(mzd2));

	GeoX = struct('mr',mxr,'ml',mxl, 'pr', xpr, 'pl', xpl);
	GeoY = struct('mf',myf,'mb',myb, 'pf', ypf, 'pb', ypb);
	GeoZ = struct('mu',mzu,'md',mzd, 'pu', zpu, 'pd', zpd);
	Geo = struct('x',GeoX,'y',GeoY,'z',GeoZ);

	% update distance map
	for i = 1:100
	
		k1 = ReInitialStep(obj, Distance,d0,Dx,Dy,Dz, Geo);
		D1 = Distance-k1;
		k2 = ReInitialStep(obj, D1,d0,Dx,Dy,Dz, Geo);
		Distance = (Distance+D1-k2)/2;

	end

	obj.F = Distance;
	disp('reinitialized ...')

end

function [Step] = ReInitialStep(obj, d, d0, Dx, Dy, Dz, Geo)

	% d 	 		current distance function
	% d0			initial distance function
	% Dx,Dy,Dz		grid spacing in x,y,z direction
	% Step 			change in Distance to be applied	

	[Ny, Nx, Nz] = size(d0);
	Epsilon = 10^(-10);

	Step = zeros(Ny,Nx,Nz); % initialize step
	deltat = min(min(Dx,Dy),Dz)*ones(Ny,Nx,Nz); % initialize the adaptive time step as minum of grid spacing
	xp = Dx*ones(Ny,Nx,Nz); % grid spacing in x direction
	yp = Dy*ones(Ny,Nx,Nz); % grid spacing in y direction
	zp = Dz*ones(Ny,Nx,Nz); % grid spacing in z direction

	% initialize derivative matrices

		xR = zeros(Ny,Nx,Nz); xL = xR; yR = xR; yL = xR; zR = xR; zL = xR;

	% ENO derivatives

%		dl = circshift(d, [0 1 0]);  dl(:,1,:)   = d(:,1,:);	% x:left
%		dr = circshift(d, [0 -1 0]); dr(:,end,:) = d(:,end,:);	% x:right
%		df = circshift(d, [-1 0 0]); df(end,:,:) = d(end,:,:);	% y:front
%		db = circshift(d, [1 0 0]);  db(1,:,:)   = d(1,:,:);	% y:back
%		du = circshift(d, [0 0 -1]); du(:,:,end) = d(:,:,end);	% z:up
%		dd = circshift(d, [0 0 1]);  dd(:,:,1)   = d(:,:,1);	% z:down

		dl = d(obj.GD3.oxo); dl(:,1,:)   = d(:,1,:);	% x:left
		dr = d(obj.GD3.oXo); dr(:,end,:) = d(:,end,:);	% x:right
		df = d(obj.GD3.Yoo); df(end,:,:) = d(end,:,:);	% y:front
		db = d(obj.GD3.yoo); db(1,:,:)   = d(1,:,:);	% y:back
		du = d(obj.GD3.ooZ); du(:,:,end) = d(:,:,end);	% z:up
		dd = d(obj.GD3.ooz); dd(:,:,1)   = d(:,:,1);	% z:down

	% Newton's second difference in x and y direction

			D2x = (dl+dr-2*d)/Dx^2/2;
			D2y = (df+db-2*d)/Dy^2/2;
			D2z = (du+dd-2*d)/Dz^2/2;

%			D2xl = circshift(D2x, [0 1 0]);  D2xl(:,1,:)   = D2x(:,1,:);
%			D2xr = circshift(D2x, [0 -1 0]); D2xr(:,end,:) = D2x(:,end,:);
%			D2yf = circshift(D2y, [-1 0 0]); D2yf(end,:,:) = D2y(end,:,:);
%			D2yb = circshift(D2y, [1 0 0]);  D2yb(1,:,:)   = D2y(1,:,:);
%			D2zu = circshift(D2z, [0 0 -1]); D2zu(:,:,end) = D2z(:,:,end);
%			D2zd = circshift(D2z, [0 0 1]);  D2zd(:,:,1)   = D2z(:,:,1);

			D2xl = D2x(obj.GD3.oxo); D2xl(:,1,:)   = D2x(:,1,:);
			D2xr = D2x(obj.GD3.oXo); D2xr(:,end,:) = D2x(:,end,:);
			D2yf = D2y(obj.GD3.Yoo); D2yf(end,:,:) = D2y(end,:,:);
			D2yb = D2y(obj.GD3.yoo); D2yb(1,:,:)   = D2y(1,:,:);
			D2zu = D2z(obj.GD3.ooZ); D2zu(:,:,end) = D2z(:,:,end);
			D2zd = D2z(obj.GD3.ooz); D2zd(:,:,1)   = D2z(:,:,1);

		% rightward and leftward x derivative

			mmodxr = MinMod(D2x, D2xr); xR = (dr-d)/Dx - Dx * mmodxr;
			mmodxl = MinMod(D2x, D2xl); xL = (d-dl)/Dx + Dx * mmodxl;

			mxr = Geo.x.mr;
			mxl = Geo.x.ml;

			xp = Geo.x.pr;
			xR(mxr) = -d(mxr)./xp(mxr) - xp(mxr).*mmodxr(mxr);

			xp(mxl) = Geo.x.pl(mxl);
			xL(mxl) = d(mxl)./xp(mxl) + xp(mxl).*mmodxl(mxl);

		% forward and backward y derivative

			mmodyf = MinMod(D2y, D2yf); yF = (df-d)/Dy - Dy * mmodyf;
			mmodyb = MinMod(D2y, D2yb); yB = (d-db)/Dy + Dy * mmodyb;

			myf = Geo.y.mf;
			myb = Geo.y.mb;

			yp = Geo.y.pf;
			yF(myf) = -d(myf)./yp(myf) - yp(myf).*mmodyf(myf);

			yp(myb) = Geo.y.pb(myb);
			yB(myb) = d(myb)./yp(myb) + yp(myb).*mmodyb(myb);

		% upward and downward z derivative

			mmodzu = MinMod(D2z, D2zu); zU = (du-d)/Dz - Dz * mmodzu;
			mmodzd = MinMod(D2z, D2zd); zD = (d-dd)/Dz + Dz * mmodzd;

			mzu = Geo.z.mu;
			mzd = Geo.z.md;

			zp = Geo.z.pu;
			zU(mzu) = -d(mzu)./zp(mzu) - zp(mzu).*mmodzu(mzu);

			zp(mzd) = Geo.z.pd(mzd);
			zD(mzd) = d(mzd)./zp(mzd) + zp(mzd).*mmodzd(mzd);

		% adpative time step

			deltat = min(deltat, abs(xp));
			deltat = min(deltat, abs(yp));
			deltat = min(deltat, abs(zp));
			deltat = 0.3 * deltat;

		% calculate the time step

			Pos = d0>=0;

			% if outside
			Step(Pos) = sqrt( max( (max(0,xL(Pos))).^2, (min(0,xR(Pos))).^2 ) ...
							+ max( (max(0,yB(Pos))).^2, (min(0,yF(Pos))).^2 ) ...
							+ max( (max(0,zD(Pos))).^2, (min(0,zU(Pos))).^2 ) );
			% if inside
			Step(~Pos) = sqrt( max( (min(0,xL(~Pos))).^2, (max(0,xR(~Pos))).^2 ) ...
							+  max( (min(0,yB(~Pos))).^2, (max(0,yF(~Pos))).^2 ) ...
							+  max( (min(0,zD(~Pos))).^2, (max(0,zU(~Pos))).^2 ) );

			Step = deltat.*sign(d0).*(Step-1);

end

function [x] = MinMod(x1, x2)
	
	a1 = abs(x1) < abs(x2);
	a2 = x1 .* x2 < 0;

	x = x2;
	x(a1) = x1(a1);
	x(a2) = 0; 

end