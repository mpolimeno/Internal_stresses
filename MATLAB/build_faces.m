function [finalposint, finalndir, finalori,Nf] = build_faces(xc, NC)
% builds the surface of the object make of cubes of side 2 with centers at location xc and number NC
% This files builds the face centers, normal directions, and orientations of a volume
% made of cubes sharing at least one face.
% Thankfully, the order does not matter.

clear('posint')

% NC = 6; % number of cubes
% xc(i,:) % the position of the centers, of the form (2*j,2*k,2*m), for j, k, m elements of Z
% xc(1,:) = [0 0 0];
% xc(2,:) = [2 0 0];
% xc(3,:) = [4 0 0];
% xc(4,:) = [2 2 0];
% xc(5,:) = [0 2 0];
% xc(6,:) = [0 0 2];


% base cube with one point per face
Nf = 6;
posintb(1,:) = [1,0,0];
posintb(2,:) = [0,1,0];
posintb(3,:) = [0,0,1];
posintb(4,:) = [0,0,-1];
posintb(5,:) = [0,-1,0];
posintb(6,:) = [-1,0,0];
ndirb(1) = 1;
orib(1) = 1;
ndirb(2) = 2;
orib(2) = 1;
ndirb(3) = 3;
orib(3) = 1;
ndirb(4) = 3;
orib(4) = -1;
ndirb(5) = 2;
orib(5) = -1;
ndirb(6) = 1;
orib(6) = -1;




posint= []; % face centers
ndir = []; % normal directions
ori = []; % normal orirentations
toremove = [];
for i=1:NC
  % add all 6 face centers
	for kk = 1:6
		posint = [posint; xc(i,:)+posintb(kk,:)];
		ndir = [ndir; ndirb(kk)];
		ori = [ori; orib(kk)];
	end;
	% check if there is an cube adjacent cube already
	for j=1:i-1
		for kk = 1:6
			if (xc(i,:) + posintb(kk,:)*2 == xc(j,:))
				toremove = [toremove; i,kk]; % remove that face from the new cube
				toremove = [toremove; j,7-kk]; % remove corresponding face from the old cube
			end;
		end;
	end;
end;

% Now remove all the useless faces
trl = length(toremove);
tokeep = [];
for j=1:6*NC
	keep = 0;
	for i=1:trl
		indx = 6*(toremove(i,1)-1)+toremove(i,2);
		if (j== indx)
			keep = 1;
		end;
	end;
	if (keep == 0)
		tokeep = [tokeep j];
	end;
end;

finalposint = posint(tokeep,:);
finalndir = ndir(tokeep);
finalori = ori(tokeep);

% that should do it, but it is completely untested.
% this is not very modular, as it if difficult to add a cube after the fact.
% it could be modified to remove a face as soon as it is found redundant

% to see the result
Nf = length(finalposint);







