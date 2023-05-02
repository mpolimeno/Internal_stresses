function [forceout,drag,torque,drag_strain,eigvec] = fractal_bi_stokes_force_strain(xc,posint,ndir,ori,Estrain,Nf);
% This function finds the force on each face of a shape moving with constant velocity in Stokes flow
% Estrain is the traceless, symmetric strain matrix, with velocity u = -E.x
% for example Estrain = [1 0 0; 0 -1 0; 0 0 0]
% we compute the flow satisfying uf = -E*x on the boundary and 0 at infinity
% and then consider u= E*x + u
% posint is the list of th ecenter location of each square face
% ndir is the direction of the normal of each face (1,2 or 3)
% ori is the orientation of each square face -1 or 1
% Nf is the number of faces
% the force found can be used to evaluate the velocity elsewhere using
% [velout] = fractal_bi_stokes_vel(posint,ndir,ori,Uvec,mypoint,Nf,forceout);

mu = 1; % for now
LHS = zeros(6*3,6*3);
RHStemp = zeros(6*3,6*3);

uall = [];
xhom = []; % this will be the homogeneous solution, where all the forces are parallel to the normals.
            % We do not want it;
sz = size(xc);
if (sz(1)==1)
	cm = xc;
else
	cm = mean(xc); % position of the center of mass
end;

for kk = 1:Nf

	v2 = (posint(kk,:) - cm)';
	Uvec = -Estrain*v2; % velocity at the center of each face under strain E

  uall = [uall ; Uvec ];
	if (ndir(kk) == 1)
		xhom = [xhom; ori(kk)*[1;0;0]];
	end;
	if (ndir(kk) == 2)
		xhom = [xhom; ori(kk)*[0;1;0]];
	end;
	if (ndir(kk) == 3)
		xhom = [xhom; ori(kk)*[0;0;1]];
	end;
end;


% build the coefficients of force by integrating on the surface
for kk = 1:Nf
	pos0 = posint(kk,:); % where we evaluate the velocity
  for ff=1:Nf
		myposint = posint(ff,:); % where we integrate
		myndir = ndir(ff); % normal of the integrated surface
		myori = ori(ff); % orientation of the normal
		[constij,xxij] = single_layer(pos0,myposint,myndir);
		SL = (constij + xxij)*(-1/(8*pi*mu)); % maybe with mu?

		% DL = double_layer(pos0,myposint,myndir,myori)*(1/(4*pi)); % maybe with mu?
 		% RHStemp(3*(ff-1)+1:3*ff,3*(kk-1)+1:3*kk) = DL;

		% LHS(3*(ff-1)+1:3*ff,3*(kk-1)+1:3*kk) = SL;
		LHS(3*(kk-1)+1:3*kk,3*(ff-1)+1:3*ff) = SL;

	end;
end;



% Expand the system to force a solution orthogonal to the solution parallel to all normals
LHS = [LHS; xhom'];
RHS = [uall;0]; % add the zero to force solution orthogonal to xhom
force = LHS\RHS;
% r = rank(LHS,eps);
% if r==size(LHS,2)
%     force = LHS\RHS;
% else
%     force = pinv(LHS)*RHS;
% end
% Now, you should be able to evaluate the velocity anywhere
forceout = [];
drag = [0;0;0];
torque = [0;0;0];
% get eigenvectors
[eigvec eigval] = eig(Estrain);

drag_strain(1:3) = 0;
% Add back the strain from the term E*x
for ff=1:Nf
	v1 = (posint(ff,:) - cm)';
	if (ndir(ff) == 1)
		nor1 = [1;0;0]*ori(ff);
	end;
	if (ndir(ff) == 2)
		nor1 = [0;1;0]*ori(ff);
	end;
	if (ndir(ff) == 3)
		nor1 = [0;0;1]*ori(ff);
	end;
	v2 = (force(3*(ff-1)+1:3*ff) + 2*Estrain*nor1); % corrected local stress
	forceout = [forceout v2]; % total stress
	drag = drag + 4*v2;
	torque = torque + 4*(cross(v1,v2));
	for kk= 1:3
		drag_strain(kk) = drag_strain(kk) + 4*abs(v2'*eigvec(:,kk));
	end;
end;

% Remove any net drag
for kk= 1:3
	drag_strain(kk) = (drag_strain(kk) - abs(drag'*eigvec(:,kk)))/2;
end;








