% Below is a sample of how to use the codes in this folder
% We compare the effects of the resolution for either the dumbell or for
% DLA aggregates

% method:
% 1 = DLA
% 2 = dumbell
method = 1;


if method==1
    % First create an aggregate via DLA
    NC = 20; % number of cubes
    % NOTE in DLA 3D, rng is seeded at SEED
    SEED = 1;
    [xc] = DLA_3D(NC,SEED);
else
    NC_dumb = 4;
    NC_bell = 54;
    pos_dumb = build_dumb(NC_dumb,1);
    [pos_bell_l,pos_bell_r] = build_bell(NC_bell,1);
    xc = [pos_bell_l;pos_dumb;pos_bell_r];
    NC = NC_dumb + NC_bell;
end

%resolution of agggregate
%res = 1 -> low res
%res = 2 -> high res
res = 1;
if res==1
    % for our kinds of flows, we need the aggregate to be centered at the
    % % origin for this formulation to hold
    if size(xc,1)==1
        cm = xc; % for one cube only
    else
        cm = mean(xc); % for multiple cubes
    end

    % Now compute where are the faces and what are their normals and orientations
    [finalposint, finalndir, finalori,Nf] = build_faces(xc, NC);
    % Nf is the number of faces
    % finalposint is the location of the center of each face
    % finalndir is the direction of the normal
    % finalori is the orientation of each face
    % Nf is the number of faces

end


if res==2
    xc_low_res = xc;
    % higher resolution dumbell
    [xc,NC] = refine_aggregate(xc_low_res, size(xc_low_res,1));
    % for our kinds of flows, we need the aggregate to be centered at the
    % % origin for this formulation to hold
    if size(xc,1)==1
        cm = xc; % for one cube only
    else
        cm = mean(xc); % for multiple cubes
    end
    
    % Now compute where are the faces and what are their normals and orientations
    [finalposint, finalndir, finalori,Nf] = build_faces(xc,NC);
    % Nf is the number of faces
    % finalposint is the location of the center of each face
    % finalndir is the direction of the normal
    % finalori is the orientation of each face
    % Nf is the number of faces
end


%%%%%%%%%%%%%%%% SELECT FLOW %%%%%%%%%%%%%%%%%%%

% all 4 of these flows recover the correct velocities
% flow = 1 -> translation
% flow = 2 -> rotation
% flow = 3 -> translation + rotation % case 3 needs to be tested more
% flow = 4 -> extensional flow % case 4 needs more testing too
% flow = 5 -> try to break dumbell in the middle
flow = 4;

if flow~=1 && flow~=2 && flow~=3 && flow~=4 && flow~=5
    msg = "Flow must be either 1, 2, 3, 4 or 5";
    error(msg);
end
% define matrix for extensional flow
% from here it will also get passed to Double Layer solver
% if the flow is of other type, pass empty vector
if flow==4
   M = [-1,0,0;0,1,0;0,0,0];
else
    M = [];
end

% based on different flows, select appropriate codes to get drag and torque
% to be fed to the new code with double layer
if flow==1
    msg = "Translational flow selected";
    disp(msg);
    
    U_infty = [0;0;1]; % translation

    % forceout is stress on each external face    
    [forceout,drag,torque] = fractal_bi_stokes_force(xc,finalposint,finalndir,finalori,U_infty,Nf);
elseif flow==2
    msg = "Rotational flow selected";
    disp(msg);
    
    U_infty = [0;0;1]; % rotation

    [forceout,drag,torque] = fractal_bi_stokes_force_rot(xc,finalposint,finalndir,finalori,U_infty,Nf);
    [Amat,bmat,allforces] = inner_stresses(xc,NC,forceout,finalposint,finalndir,finalori,U_infty,Nf,drag);
elseif flow==3
    msg = "Translational + Rotational flow selected";
    disp(msg);
    
    Uvec = [0;0;1];
    Rot = [0;0;1];
    
    [forceout,drag,torque] = fractal_bi_stokes_force_trans_and_rot(xc,finalposint,finalndir,finalori,Uvec,Rot,Nf);
elseif flow==4
    msg = "Extensional flow selected";
    disp(msg);
    % extensional flow
    
    % drag_strain should not be needed to solve linear system
    [forceout,drag,torque,drag_strain,eigvec] = fractal_bi_stokes_force_strain(xc,finalposint,finalndir,finalori,M,Nf);
else
    msg_now = "ERROR: FLOW NOT READY YET";
    error(msg_now);
end

% fixing drag and torque
drag_in = drag;
torque_in = torque; 

tic
% Compute External Stresses and Velocities for solid body motion
[LHS,sol,stress_outer,U_vec,Omega_vec] = ComputeStressesAndSolidBodyMotion(xc,finalposint,finalndir,finalori,drag_in,torque_in,Nf,M,flow);

% Use external stresses to compute internal stresses
[internal_and_external_stresses,internal_stresses,indices_of_internal_stresses,internal_faces_and_cubes_index_array,internal_faces_and_cubes_index_array_no_double_counting] = ComputeInternalStresses(xc,NC,stress_outer,U_vec,Omega_vec,drag_in);
toc       

T = table(internal_stresses')
filename = sprintf('Internal_Stresses_%i.txt',SEED)
%writetable(T,'res.txt','Delimiter','\t','WriteRowNames',false)
writetable(T,filename,'Delimiter','\t','WriteRowNames',false)
