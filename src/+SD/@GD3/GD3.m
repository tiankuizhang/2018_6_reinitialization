classdef GD3 < handle
	
	%Grid3 : meshgrid in 3D

	properties (SetAccess = immutable)
		X % meshgrid X
		Y % meshgrid Y
		Z % meshgrid Z
		Dx % spacing in x direction
		Dy % spacing in y direction
		Dz % spacing in z direction
		Ds % mean spacing
		EltVol % volume of one element
		ncols % number of grid points in x direction
		mrows % number of grid points in y direction
		lshts % number of grid points in z direction
		Size % size array
		NumElt % total number of grid points 
		xmin % min value of coordinate x 
		xmax % max value of coordinate x
		ymin % min value of coordinate y
		ymax % max value of coordinate y
		zmin % min value of coordinate z
		zmax % max value of coordinate z
		BOX % bounding box of grid domain
		PX	% permutation of X for griddedInterpolant
		PY	% permutation of Y for griddedInterpolant
		PZ	% permutation of Z for griddedInterpolant
		
		% index matrix and the shifted index matrices
		ooo % not shifted
		
		oXo % shifted in the -x direction thus represent indices to the right  
		oxo % shifted in the +x direction thus represent indices to the left
		Yoo 
		yoo
		ooZ
		ooz

		soXo % smae with above definition but without periodic boundary condtion
		soxo
		sYoo
		syoo
		sooZ
		sooz
		
		YXo
		yXo
		Yxo
		yxo
		
		YoZ
		Yoz
		yoZ
		yoz
		
		oXZ
		oxZ
		oXz
		oxz
		
		YXZ
		YXz
		YxZ
		Yxz
		yXZ
		yXz
		yxZ
		yxz

		% first/sencond order derivative operators (including cross derivatives)
		% using centeral difference with periodic boundary condition
		Lx
		Ly
		Lz
		Lxx
		Lyy
		Lzz
		Lxy
		Lyz
		Lzx

		% identity matrix
		Idt

	end

	methods

		function obj = GD3(Xm, Ym, Zm)
			obj.X = Xm;
			obj.Y = Ym;
			obj.Z = Zm;
			obj.Dx = Xm(1,2,1) - Xm(1,1,1);
			obj.Dy = Ym(2,1,1) - Ym(1,1,1);
			obj.Dz = Zm(1,1,2) - Zm(1,1,1);
			obj.Ds = (obj.Dx + obj.Dy + obj.Dz) / 3;
			obj.EltVol = obj.Dx * obj.Dy * obj.Dz;
			obj.Size = size(Xm);
			obj.mrows = obj.Size(1);
			obj.ncols = obj.Size(2);
			obj.lshts = obj.Size(3);
			obj.NumElt = prod(obj.Size);
			obj.xmin = Xm(1,1,1);
			obj.ymin = Ym(1,1,1);
			obj.zmin = Zm(1,1,1);
			obj.xmax = Xm(1,end,1);
			obj.ymax = Ym(end,1,1);
			obj.zmax = Zm(1,1,end);
			obj.BOX = [obj.xmin obj.xmax obj.ymin obj.ymax obj.zmin obj.zmax];
			obj.PX = permute(Xm, [2 1 3]);
			obj.PY = permute(Ym, [2 1 3]);
			obj.PZ = permute(Zm, [2 1 3]);

			% the largest integer uint32 can represent is 2^32~4e9~(1.6e3)^3
			% the largest integer uint16 can represent is 2^16~6e4~(4e2)^40^3
			% thus as long as NumElt<2^32 we can use uint32 
			obj.ooo = reshape(1:obj.NumElt,obj.Size);

			obj.oXo = circshift(obj.ooo, [	0	-1	0	]);
			obj.oxo = circshift(obj.ooo, [	0	1	0	]);
			obj.Yoo = circshift(obj.ooo, [	-1	0	0	]);
			obj.yoo = circshift(obj.ooo, [	1	0	0	]);
			obj.ooZ = circshift(obj.ooo, [	0	0	-1	]);
			obj.ooz = circshift(obj.ooo, [	0	0	1	]);

			obj.soXo = obj.oXo; obj.soXo(:,end,:) = obj.ooo(:,end,:);
			obj.soxo = obj.oxo; obj.soxo(:,1  ,:) = obj.ooo(:,1  ,:);
			obj.sYoo = obj.Yoo; obj.sYoo(end,:,:) = obj.ooo(end,:,:);
			obj.syoo = obj.yoo; obj.syoo(1  ,:,:) = obj.ooo(1  ,:,:);
			obj.sooZ = obj.ooZ; obj.sooZ(:,:,end) = obj.ooo(:,:,end);
			obj.sooz = obj.ooz; obj.sooz(:,:,1  ) = obj.ooo(:,:,1   );

			obj.YXo = circshift(obj.ooo, [	-1	-1	0	]);
			obj.yXo = circshift(obj.ooo, [	1	-1	0	]); 
			obj.Yxo = circshift(obj.ooo, [	-1	1	0	]);
			obj.yxo = circshift(obj.ooo, [	1	1	0	]); 

			obj.YoZ = circshift(obj.ooo, [	-1	0	-1	]);
			obj.Yoz = circshift(obj.ooo, [	-1	0	1	]); 
			obj.yoZ = circshift(obj.ooo, [	1	0	-1	]);
			obj.yoz = circshift(obj.ooo, [	1	0	1	]); 

			obj.oXZ = circshift(obj.ooo, [	0	-1	-1	]);
			obj.oxZ = circshift(obj.ooo, [	0	1	-1	]); 
			obj.oXz = circshift(obj.ooo, [	0	-1	1	]);
			obj.oxz = circshift(obj.ooo, [	0	1	1	]); 

			obj.YXZ = circshift(obj.ooo, [	-1	-1	-1	]);
			obj.YXz = circshift(obj.ooo, [	-1	-1	1	]); 
			obj.YxZ = circshift(obj.ooo, [	-1	1	-1	]);
			obj.Yxz = circshift(obj.ooo, [	-1	1	1	]); 
			obj.yXZ = circshift(obj.ooo, [	1	-1	-1	]);
			obj.yXz = circshift(obj.ooo, [	1	-1	1	]); 
			obj.yxZ = circshift(obj.ooo, [	1	1	-1	]);
			obj.yxz = circshift(obj.ooo, [	1	1	1	]); 

			obj.Lx  =	( sparse(obj.ooo(:), obj.oXo(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.oxo(:),	-1,	obj.NumElt, obj.NumElt) ) / (2*obj.Dx);

			obj.Ly  =	( sparse(obj.ooo(:), obj.Yoo(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.yoo(:),	-1,	obj.NumElt, obj.NumElt) ) / (2*obj.Dy);

			obj.Lz  = 	( sparse(obj.ooo(:), obj.ooZ(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.ooz(:),	-1,	obj.NumElt, obj.NumElt) ) / (2*obj.Dz); 

			obj.Lxx =	( sparse(obj.ooo(:), obj.oXo(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.ooo(:),	-2,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.oxo(:),	1,	obj.NumElt, obj.NumElt) ) / (obj.Dx.^2);

			obj.Lyy =	( sparse(obj.ooo(:), obj.Yoo(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.ooo(:),	-2,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.yoo(:),	1,	obj.NumElt, obj.NumElt) ) / (obj.Dy.^2);

			obj.Lzz = 	( sparse(obj.ooo(:), obj.ooZ(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.ooo(:),	-2,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.ooz(:),	1,	obj.NumElt, obj.NumElt) ) / (obj.Dz.^2);

			obj.Lxy = 	( sparse(obj.ooo(:), obj.YXo(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.yxo(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.Yxo(:),	-1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.yXo(:),	-1,	obj.NumElt, obj.NumElt) ) / (4*obj.Ds.^2); 

			obj.Lyz = 	( sparse(obj.ooo(:), obj.YoZ(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.yoz(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.Yoz(:),	-1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.yoZ(:),	-1,	obj.NumElt, obj.NumElt) ) / (4*obj.Ds.^2); 

			obj.Lzx = 	( sparse(obj.ooo(:), obj.oXZ(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.oxz(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.oXz(:),	-1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.oxZ(:),	-1,	obj.NumElt, obj.NumElt) ) / (4*obj.Ds.^2); 

			obj.Idt = sparse(obj.ooo(:), obj.ooo(:), 1, obj.NumElt, obj.NumElt);
		end

		% create a sparse diagonal matrix out of a field
		% look up the function spdiags and replace it with the built-in function
		function val = SparseDiag(obj, Field)
			val = sparse(obj.ooo(:), obj.ooo(:), reshape(Field,[],1), obj.NumElt, obj.NumElt);
		end

		% various derivatives of a field
		function val = Fx(obj, Field)
			val = (Field(obj.oXo) - Field(obj.oxo)) / (2*obj.Dx);
		end

		function val = Fy(obj, Field)
			val = (Field(obj.Yoo) - Field(obj.yoo)) / (2*obj.Dy);
		end

		function val = Fz(obj, Field)
			val = (Field(obj.ooZ) - Field(obj.ooz)) / (2*obj.Dz);
		end

		function val = Fxx(obj, Field)
			val = (Field(obj.oXo) - 2*Field + Field(obj.oxo)) / (obj.Dx.^2);
		end

		function val = Fyy(obj, Field)
			val = (Field(obj.Yoo) - 2*Field + Field(obj.yoo)) / (obj.Dy.^2);
		end

		function val = Fzz(obj, Field)
			val = (Field(obj.ooZ) - 2*Field + Field(obj.ooz)) / (obj.Dz.^2);
		end

		function val = Fxy(obj, Field)
			val = (Field(obj.YXo) + Field(obj.yxo) - Field(obj.Yxo) - Field(obj.yXo)) / (4*obj.Ds.^2);
		end

		function val = Fyz(obj, Field)
			val = (Field(obj.YoZ) + Field(obj.yoz) - Field(obj.Yoz) - Field(obj.yoZ)) / (4*obj.Ds.^2);
		end

		function val = Fzx(obj, Field)
			val = (Field(obj.oXZ) + Field(obj.oxz) - Field(obj.oXz) - Field(obj.oxZ)) / (4*obj.Ds.^2);
		end

		function val = Laplacian(obj, Field)
			val = obj.Fxx(Field) + obj.Fyy(Field) + obj.Fzz(Field);
		end


	end

end