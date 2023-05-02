function [xc] = DLA_3D(NC,seeding)
% This function generates a 3D DLA cluster
% NC is the number of cubes to place
% THIS MAY HAVE A BUG AS IT SOMETIMES ADDS A CUBE WHERE A CUBE IS ALREADY PRESENT

 rng(seeding); % seed the randomness

if ((NC < 1) || (round(NC) ~= NC))
	error('Are you kidding? I need a positive integer as an input');
end;

xc = zeros(NC,3);
xc(1,:) = [0 0 0]; % start at the origin %seems repetitive, it is initialized to zero

for i=2:NC
	%cm = mean(xc(1:i-1)); % find center of mass
    cm = mean(xc(1:i-1,:)); % find center of mass
    % shouldn't it be cm = mean(xc(1:i-1,:))?
	for j=1:i-1
		radj(j) = norm(xc(j,:)-cm); % distance to center
	end;
	rad = max(radj); % size of the cluster
	myrad = rad+10; % radius of starting point

	success = 0;
	while (success == 0)
		% this phi (polar angle) does not guarantee uniform distribution
        % see here: https://mathworld.wolfram.com/SpherePointPicking.html
        %phi = rand*pi;
        
        phi = acos(2*rand - 1); % correct phi
		theta = rand*2*pi;
		posrw = cm + [(myrad)*cos(theta)*sin(phi),myrad*sin(theta)*sin(phi),myrad*cos(phi)];
		posrw = 2*round(posrw/2); % true starting point

		% random walk until you go too far or find the cluster
		wearedone = 0;
		while (wearedone == 0)
			ra = rand; % to pick a direction of motion
			% now take a step
			if ((0 <= ra) && (ra < 0.1666666667))
				posrw = posrw + [2 0 0];
			end;
			if ((0.1666666667 <= ra) && (ra < 0.333333333))
				posrw = posrw + [-2 0 0];
			end;
			if ((0.333333333 <= ra) && (ra < 0.5))
				posrw = posrw + [0 2 0];
			end;
			if ((0.5 <= ra) && (ra < 0.6666666667))
				posrw = posrw + [0 -2 0];
			end;
			if ((0.6666666667 <= ra) && (ra < 0.833333333))
				posrw = posrw + [0 0 2];
			end;
			if ((0.833333333 <= ra) && (ra < 1))
				posrw = posrw + [0 0 -2];
			end;
			% check how far we are
			newrad = norm(posrw-cm);
			if (newrad > myrad+10)
				wearedone = 1;
				success = 0;
			end;
		% check who our neighbours are
			if (wearedone == 0)
				for j=1:i-1
					if (norm(posrw - xc(j,:)) == 2)
						wearedone = 1;
						success = 1;
					end;
				end;
			end;
		end;
	end;
	xc(i,:) = posrw;
end;

% NC3 = 186;
% xc3 = xc(1:NC3,:);
% [finalposint3, finalndir3, finalori3,Nf3] = build_faces(xc3, NC3);
% [forceout3,drag3] = fractal_bi_stokes_force(finalposint3,finalndir3,finalori3,Uvec,Nf3);


