function [forceout,drag,torque] = fractal_bi_stokes_force(xc,posint,ndir,ori,Uvec,Nf);
%function forceout = fractal_bi_stokes_force(xc,posint,ndir,ori,Uvec,Nf);
% This function finds the force on each face of a shape moving with constant velocity in Stokes flow
% Uvec is the direction of motion of the object, should be column vector
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
            % We do not want it
for kk = 1:Nf
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
%disp(size(xhom))
% build the coefficients of force by integrating on the surface
for kk = 1:Nf
	pos0 = posint(kk,:); % where we evaluate the velocity
  for ff=1:Nf
		myposint = posint(ff,:); % where we integrate
		myndir = ndir(ff); % normal of the integrated surface
		myori = ori(ff); % orientation of the normal
		[constij,xxij] = single_layer(pos0,myposint,myndir); % this builds the matrices

        %[intij] = double_layer(pos0,myposint,myndir,myori); % just prrint it to see it
		SL = (constij + xxij)*(-1/(8*pi*mu)); % maybe with mu?
   
		% DL = double_layer(pos0,myposint,myndir,myori)*(1/(4*pi)); % maybe with mu?
 		% RHStemp(3*(ff-1)+1:3*ff,3*(kk-1)+1:3*kk) = DL;

		% LHS(3*(ff-1)+1:3*ff,3*(kk-1)+1:3*kk) = SL;
		LHS(3*(kk-1)+1:3*kk,3*(ff-1)+1:3*ff) = SL;

	end;
end;

% eig(RHStemp*8*pi)
% RHStemp*ones(Nf*3,1)
% pause
%disp(size(LHS))
% Expand the system to force a solution orthogonal to the solution parallel to all normals
LHS = [LHS; xhom'];
%disp(size(LHS))
RHS = [uall;0]; % add the zero to force solution orthogonal to xhom
%disp(size(RHS))
uall;
debugout = LHS(1:18,1:18);

% solve system

force = LHS\RHS; %Matteo: These are the external stresses

%disp(force)
% Now, you should be able to evaluate the velocity anywhere
forceout = [];
drag = [0;0;0];
torque = [0;0;0];
cm = mean(xc);
for ff=1:Nf
	forceout = [forceout force(3*(ff-1)+1:3*ff)];
 	drag = drag + 4*force(3*(ff-1)+1:3*ff);
 	v1 = (posint(ff,:) - cm)';
 	v2 = force(3*(ff-1)+1:3*ff);
 	torque = torque + 4*(cross(v1,v2));
end;
