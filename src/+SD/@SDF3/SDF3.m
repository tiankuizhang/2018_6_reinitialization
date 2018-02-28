classdef SDF3 < handle
	
	%SDF3 : signed distance function in 3D

	properties (SetAccess = immutable)
		GD3 % SD.GD3 object
		Sme_Thk % smeared thickness for approximating Dirac delta function and Heaviside function
		En_Volume % enclosed volume of the vesicle
		Suf_Area % surface area of the vesicle
		Red_Vol % reduced volume of the vesicle
	end

	properties 
		F % values of the signed distance function

		% components of gradients
		Fx % 1st derivative along x-axis with central difference
		Fy % 1st derivative along y-axis with central difference
		Fz % 1st derivative along z-axis with central difference
		% components of Hessian matrix
		Fxx % 2nd derivative along x-axis with central difference
		Fyy % 2nd derivative along y-axis with central difference
		Fzz % 2nd derivative along z-axis with central difference
		Fxy % cross derivatives
		Fyz % cross derivatives
		Fzx % cross derivatives

		% normal vector
		epsilon = 1e-4; % see 2003_Smereka
		Fg % magnitude of gradient with soomthing
		Fg_1 % magnitude of gradient without soomthing
		Nx % x component of normal vectors
		Ny % y component of normal vectors
		Nz % z component of normal vectors

		SC % sum of principal curvatures. unit sphere has SC equal to 2
		GC % Gaussian curvature. unit sphere has GC equal to 1
	
		Dirac_Delta % Dirac_Delta function for calculating surface integral
		Heaviside % Heaviside function for calculating volume integral

		LSL % linear part of the surface Laplacian operator
		LCF % linear part of the curvature force operator
		
		
		CF % curvature force
		NCF % nonlinear part of the curvature force 

		LND % normal derivative operator
		CF_Op % curvature force operator

	end


	methods

		function obj = SDF3(Xm, Ym, Zm, Val)
			obj.GD3 = SD.GD3(Xm, Ym, Zm);
			obj.Sme_Thk = obj.GD3.Dx * 1.5;
			%obj.reinitialization(Val)
			obj.En_Volume = obj.VolumeIntegral(1);
			obj.Suf_Area = obj.SurfaceIntegral(1);
			obj.Red_Vol = (3*obj.En_Volume/4/pi) * (4*pi/obj.Suf_Area)^(3/2);
		end

		function obj = set.F(obj, val)
			obj.F = val;

			obj.Fx = obj.GD3.Fx(val);
			obj.Fy = obj.GD3.Fy(val);
			obj.Fz = obj.GD3.Fz(val);
			obj.Fxx = obj.GD3.Fxx(val);
			obj.Fyy = obj.GD3.Fyy(val);
			obj.Fzz = obj.GD3.Fzz(val);
			obj.Fxy = obj.GD3.Fxy(val);
			obj.Fyz = obj.GD3.Fyz(val);
			obj.Fzx = obj.GD3.Fzx(val);

			obj.Fg_1 = sqrt(obj.Fx.^2 + obj.Fy.^2 + obj.Fz.^2 );
			obj.Fg = sqrt(obj.Fg_1.^2 + obj.epsilon);
			obj.Nx = obj.Fx ./ obj.Fg;
			obj.Ny = obj.Fy ./ obj.Fg;
			obj.Nz = obj.Fz ./ obj.Fg;

		%	obj.LSL = obj.GD3.Lxx + obj.GD3.Lyy + obj.GD3.Lzz ...
		%			- ( obj.GD3.SparseDiag(obj.Nx .* obj.Nx) * obj.GD3.Lxx + ...
		%				obj.GD3.SparseDiag(obj.Ny .* obj.Ny) * obj.GD3.Lyy + ...
		%				obj.GD3.SparseDiag(obj.Nz .* obj.Nz) * obj.GD3.Lzz + ...
		%				obj.GD3.SparseDiag(obj.Nx .* obj.Ny) * obj.GD3.Lxy * 2 + ...
		%				obj.GD3.SparseDiag(obj.Ny .* obj.Nz) * obj.GD3.Lyz * 2 + ...
		%				obj.GD3.SparseDiag(obj.Nz .* obj.Nx) * obj.GD3.Lzx * 2  );

		%	obj.LND = obj.GD3.SparseDiag(obj.Nx) * obj.GD3.Lx + ...
		%			  obj.GD3.SparseDiag(obj.Ny) * obj.GD3.Ly + ...
		%			  obj.GD3.SparseDiag(obj.Nz) * obj.GD3.Lz ;


			%obj.LCF = obj.GD3.SparseDiag(obj.Fg_1) * obj.LSL * obj.GD3.SparseDiag(1./obj.Fg) * obj.LSL;
						
			obj.SC = (obj.Fxx + obj.Fyy + obj.Fzz ...
				- obj.Nx .* obj.Nx .* obj.Fxx ... 
				- obj.Ny .* obj.Ny .* obj.Fyy ...
				- obj.Nz .* obj.Nz .* obj.Fzz ...
				- obj.Nx .* obj.Ny .* obj.Fxy * 2 ...
				- obj.Ny .* obj.Nz .* obj.Fyz * 2 ...
				- obj.Nz .* obj.Nx .* obj.Fzx * 2 ) ./ obj.Fg;

			% comatrix of the hessian matrix, used to calculate the Gaussian curvature
			CoHessian11 = obj.Fyy .* obj.Fzz - obj.Fyz.^2;
			CoHessian22 = obj.Fzz .* obj.Fxx - obj.Fzx.^2;
			CoHessian33 = obj.Fxx .* obj.Fyy - obj.Fxy.^2;
			
			CoHessian12 = obj.Fyz .* obj.Fzx - obj.Fxy .* obj.Fzz;
			CoHessian13 = obj.Fxy .* obj.Fyz - obj.Fzx .* obj.Fyy;
			CoHessian23 = obj.Fxy .* obj.Fzx - obj.Fyz .* obj.Fxx;
			
			CoHessian21 = CoHessian12;
			CoHessian31 = CoHessian13;
			CoHessian32 = CoHessian23;

			obj.GC = (obj.Fx .* (obj.Fx .* CoHessian11 + obj.Fy .* CoHessian12 + obj.Fz .* CoHessian13) ...
					+ obj.Fy .* (obj.Fx .* CoHessian21 + obj.Fy .* CoHessian22 + obj.Fz .* CoHessian23) ...
					+ obj.Fz .* (obj.Fx .* CoHessian31 + obj.Fy .* CoHessian32 + obj.Fz .* CoHessian33))./obj.Fg.^4;

			% smeared Dirac delta function & Heaviside function
			In_vsc = val < -obj.Sme_Thk;
			Ot_vsc = val > obj.Sme_Thk;
			Ed_vsc = ~In_vsc & ~Ot_vsc;

			obj.Dirac_Delta = zeros(obj.GD3.Size);
			obj.Dirac_Delta(Ed_vsc) = obj.Fg_1(Ed_vsc) .* (1 + cos(pi*val(Ed_vsc)/obj.Sme_Thk)) / obj.Sme_Thk / 2;

			obj.Heaviside = zeros(obj.GD3.Size);
			obj.Heaviside(Ot_vsc) = 1;
			obj.Heaviside(Ed_vsc) = (1 + val(Ed_vsc)/obj.Sme_Thk + sin(pi*val(Ed_vsc)/obj.Sme_Thk)/pi ) / 2;

			%obj.NCF = obj.Fg_1 .* obj.SC .* obj.ND(obj.SC);
			%obj.CF = obj.SL(obj.SC);

			
			
			%obj.TotalCurvature = reshape(obj.SparseDiag(1./obj.Fg) * obj.LSL * obj.F(:), obj.GD3.Size);

		%	obj.CF_Op = obj.GD3.SparseDiag(obj.Fg_1) * ...
		%			  ( obj.LSL - obj.GD3.SparseDiag(obj.SC) * obj.LND ) * ...
		%			  obj.GD3.SparseDiag(1./obj.Fg) * obj.LSL ;
		end

	end

	methods
		% derivative in the normal direction
		function val = ND(obj, Field)
			val = obj.Nx .* obj.GD3.Fx(Field) + obj.Ny .* obj.GD3.Fy(Field) + obj.Nz .* obj.GD3.Fz(Field);
		end

		% surface Laplacian of a Field
		function val = SL(obj, Field)
			val = obj.GD3.Laplacian(Field) - obj.SC .* obj.ND(Field) ...
				- ( obj.Nx .* obj.Nx .* obj.GD3.Fxx(Field) + ...
					obj.Ny .* obj.Ny .* obj.GD3.Fyy(Field) + ...
					obj.Nz .* obj.Nz .* obj.GD3.Fzz(Field) + ...
					obj.Nx .* obj.Ny .* obj.GD3.Fxy(Field) * 2 + ...
					obj.Ny .* obj.Nz .* obj.GD3.Fyz(Field) * 2 + ...
					obj.Nz .* obj.Nx .* obj.GD3.Fzx(Field) * 2 ) ;
		end

		function val = CurvatureForce(obj)
			val = obj.SL(obj.SC);
		end

		function val = Interp(obj, Field)
			val = griddedInterpolant(obj.GD3.PX,obj.GD3.PY,obj.GD3.PZ,permute(Field,[2,1,3]),'linear' );
		end

		function val = Surf(obj, isovalue)
			val = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,obj.F, isovalue);
		end

		% surface integral of Field. Field being 1 gives surface area.
		function val = SurfaceIntegral(obj, Field)
			tmp = obj.Dirac_Delta .* Field;
			val = sum(tmp(:)) * obj.GD3.EltVol;
		end

		% volume integral of Field. Field being 1 gives volume 
		function val = VolumeIntegral(obj, Field)
			tmp = (1 - obj.Heaviside) .* Field;
			val = sum(tmp(:)) * obj.GD3.EltVol;
		end

	end

	methods
		reinitialization(obj, Distance)
	end

	methods % visualization methods

		% plot a 3D field on the val contour of the distance function
		function plotField(obj,val,Field)
			% calculate an interpolant of the Field
			Interp = griddedInterpolant(obj.GD3.PX,obj.GD3.PY,obj.GD3.PZ,permute(Field,[2,1,3]),'linear' );
			% triangle mesh of the val isosurface. TriMesh is a structure with fields "vertices" and "faces"
			TriMesh = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,obj.F,val);
			% interpolate the values of the Field onto the vertices of the triangle mesh
			SurfField = Interp(TriMesh.vertices(:,1),TriMesh.vertices(:,2),TriMesh.vertices(:,3));
			% plot surface mesh 
			patch('Vertices',TriMesh.vertices,'Faces',TriMesh.faces,'FaceVertexCData',SurfField,...
				'FaceColor','interp','EdgeColor','none')
			axis equal
			view(45,30)
			colorbar
		end

		% plot several half contours of the distance function
		function plot(obj)
			axis(obj.GD3.BOX)
			obj.plotIso(-12*obj.GD3.Dx,0.8,'Red')
			obj.plotIso(-6*obj.GD3.Dx,0.8,'Green')
			obj.plotIso(0,0.8,'Blue')
			obj.plotIso(6*obj.GD3.Dx,0.8,'Green')
			obj.plotIso(12*obj.GD3.Dx,0.8,'Red')
			daspect([1 1 1])
			view(3); 
			camlight; lighting gouraud
		end

		% plot half of the val contour of the distance function
		function plotIso(obj,val,trans,Color)
			F = obj.F;
			F(obj.GD3.Y<0) = inf;
			surf1 = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,val);
			p1 = patch(surf1);
			isonormals(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,p1)
			set(p1,'FaceColor',Color,'EdgeColor','none','FaceAlpha',trans);
		end

		% plot the val contour of the distance function
		function plotSurface(obj,val,trans,Color, time)
			F = obj.F;
			surf1 = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,val);
			p1 = patch(surf1);
			isonormals(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,p1)
			set(p1,'FaceColor',Color,'EdgeColor','k','FaceAlpha',trans);
			axis(obj.GD3.BOX)
			daspect([1 1 1])
			view(3); 
			camlight; lighting gouraud
		end
	end

end