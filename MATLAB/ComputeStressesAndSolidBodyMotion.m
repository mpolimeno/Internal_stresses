%function [stress_out,U_vec,Omega_vec] = draft_code(xc,posint,ndir,ori,drag_in,torque_in,Nf);
function [LHS,sol,stress_outer,U_vec,Omega_vec] = ComputeStressesAndSolidBodyMotion(xc,posint,ndir,ori,drag_in,torque_in,Nf,M,flow);
% This function finds the force on each face of a shape moving with constant velocity in Stokes flow
% Uvec is the direction of motion of the object, should be column vector
% posint is the list of the center location of each square face
% ndir is the direction of the normal of each face (1,2 or 3)
% ori is the orientation of each square face -1 or 1
% Nf is the number of faces

mu = 1; % for now
LHS_single = zeros(6*3,6*3);

I = [];
for ii = 1:Nf
    I = [I;eye(3)];
end

a = -1;
b = 1;
N = 10; % evaluation points for each face for the double layer
dx = (b-a)/N;
dy = dx;

aij = linspace(a+dx/2,b-dx/2,N);
[ai,bj] = meshgrid(aij,aij);

if size(xc,1)==1
    cm = xc; % for one cube only
else
    cm = mean(xc); % for multiple cubes
end

% This will contain the vector coming from the Double Layer Integral
DL_xs = [];
% % build the coefficients of force by integrating on the surface
for kk = 1:Nf
	pos0 = posint(kk,:); % where we evaluate the velocity
    
    % background velocity at x0
    if flow==1 % translational flow
        U_infty = ([0;0;1]);
    elseif flow==2 % rotation
        U_infty = (cross([0;0;1]',(pos0-cm)))';
    elseif flow==3 % translation + rotation
        U_infty = ([0;0;1]' + cross([0;0;1]',(pos0-cm)))';
    elseif flow==4
        U_infty = M*(pos0-cm)'; % extensional flow at x0
    else
        if pos0(:,2)>=0
            U_infty = (cross([0;0;1]',(pos0-cm)))';
        else
            U_infty = (cross([0;0;-1]',(pos0-cm)))';
        end  
    end
    
    DL = [0;0;0];
    for ff=1:Nf
        myposint = posint(ff,:); % where we integrate
        myndir = ndir(ff); % normal of the integrated surface
        myori = ori(ff); % orientation of the normal

        % Single Layer is solved exactly and it will be part of the LHS of
        % the linear system
        [constij,xxij] = single_layer(pos0,myposint,myndir); % this builds the matrices
        SL = (constij + xxij)*(1/(8*pi*mu));
        LHS_single(3*(kk-1)+1:3*kk,3*(ff-1)+1:3*ff) = SL;

        % Numerical Double Layer code; it will be on the RHS of the linear
        % system
        DL = DL + DoubleLayer_DL(ai,bj,myposint,myndir,myori,pos0,M,flow);

    end
    
    DL = (-(6/(8*pi))*DL*dx*dy);
    RHS_DL = DL + 0.5*U_infty;
    DL_xs = [DL_xs;RHS_DL'];
end

u_sol = [];
% this will be the homogeneous solution, where all the forces are parallel to the normals.
% We do not want it
xhom = [];

% Build the Right-Hand side
for kk = 1:Nf
    u_sol = [u_sol;DL_xs(kk,:)'];
    if (ndir(kk) == 1)
        xhom = [xhom; ori(kk)*[1;0;0]];
    elseif (ndir(kk) == 2)
        xhom = [xhom; ori(kk)*[0;1;0]];
    elseif (ndir(kk) == 3)
        xhom = [xhom; ori(kk)*[0;0;1]];
    end
end

% Build the cross-product operator for matrix on LHS
x_cross = [];
cross_op = zeros(3);
for ff=1:Nf
 	x_vec = (posint(ff,:) - cm);
    
    cross_op(1,2) = - x_vec(3);
    cross_op(2,1) =   x_vec(3);
    
    cross_op(1,3) =   x_vec(2);
    cross_op(3,1) = - x_vec(2);
    
    cross_op(2,3) = - x_vec(1);
    cross_op(3,2) =   x_vec(1);
    
    x_cross = [x_cross;cross_op];
end
% Making sure all the terms on the LHS have the right sign and appropriate
% 'orientation'
x_cross_m = -x_cross; % cross product is anticommutative
LHS = [LHS_single,I,x_cross_m];
I_t = I';
I_add = zeros(3,length(I)+6);
I_add(:,1:length(I)) = I_t;
LHS = [LHS;I_add];
x_cross_m_t = x_cross_m';
x_add = zeros(3,length(x_cross_m)+6);
x_add(:,1:length(x_cross_m)) = x_cross_m_t;
LHS = [LHS;x_add];

% % Expand the system to force a solution orthogonal to the solution parallel to all normals
x_aug = zeros(1,length(xhom)+6);
x_aug(1:length(xhom)) = xhom;
LHS = [LHS;x_aug];

% Right-Hand Side includes DL, Drag and Torque
% the factor of 4 is the area of a square face of side length 2
RHS = [u_sol;drag_in/4;torque_in/4];
RHS = [RHS;0]; % add the zero to force solution orthogonal to xhom

% Now we solve the linear system
% let us try to deal with the rank deficient issue
r = rank(LHS,eps);
if r==size(LHS,2)
    sol = LHS\RHS;
else
    sol = pinv(LHS)*RHS;
end

% Extract stresses, translational velocity and rotational velocity from
% rigid body motion
cut = 3*Nf;
stress_out = sol(1:cut);
U_vec = sol(cut+1:cut+3);
Omega_vec = sol(cut+4:end);

% store stresses in appropriate manner to be fed to the function that
% compute internal stresses
stress_outer = [];
for ff=1:Nf
	stress_outer = [stress_outer stress_out(3*(ff-1)+1:3*ff)];
end

end
