function [constij,xxij] = single_layer(pos0,posint,ndir);
% pos0 is the point where we want the velocity
% posint is the center of the of the face over which we integrate
% we want to integrate, with r = pos0 - posint
% I/r + r_i r_j/|r|^3
% the results is a matrix with indices i and j
% Pozrikidis equation 4.1.1 where G is given by 2.2.8

%Matteo: ndir is the normal to a given face?

% Matteo: square faces are mapped to (x,y,0). Integrals are then all from
% -1 to 1 for both x and y (assuming cubes of side 2)
% so all of these are built as in Appendix A of paper

% pos is r0 - r
pos = (pos0-posint);
for i=1:3
	for j=1:3

    constij(i,j) = 0;

        % Matteo: point at which we evaluate velocity is the same as the
        % face center
		if (norm(pos) == 0) % we are at the singularity,
			if (i==j)
				constij(i,j) = 8*asinh(1); %Matteo: RHS of Eq. 39 (Eunji's paper)
                % these are the diagonal terms of int(I/||x-x0||)
			else
				constij(i,j) = 0; % Matteo: non-diagonal terms
			end;
            if ((i==j) && (i~=ndir)) % Matteo: away from singulartiy
				xxij(i,j) =  4*asinh(1);
			else
				xxij(i,j) = 0;
			end;
		else % we are not at the singularity 
			% Do the term I/r first
			if (i==j)
				indx0 = mod(ndir+1,3);
				if (indx0 == 0)
					indx0 = 3;
				end;
				indy0 = mod(ndir+2,3);
				if (indy0 == 0)
					indy0 = 3;
				end;
				x0 = pos(indx0); % not ndir
				y0 = pos(indy0); % not ndir, or x0
				z0 = pos(ndir);
				z = 0; % because we took the difference already between the 2 vectors

				% first do the constant term % integrate 1/((x-x0)^2 + (y-y0)^2 + (z-z0)^2)^(1/2) dx dy
				flag = 0;
				if ((z==z0) && (abs(x0) == 1) && (abs(y0) == 1))
					flag = 1; % constant treated separately
				else
					if ((abs(x0) ~= 1) && (abs(y0) ~= 1))
						if (z==z0)
							% integrate 1/(x^2+y^2)^(1/2) dx =  log(sqrt(x^2 + y^2) + x)
							% integrate  log(sqrt((x)^2 + y^2) + x)) dy = y log(sqrt(x^2 + y^2) + x) + x log(sqrt(x^2 + y^2) + y) - y

							p1 = (1-y0)*log(sqrt((1-x0)^2 + (1-y0)^2) + (1-x0)) + (1-x0)*log(sqrt((1-x0)^2 + (1-y0)^2) + (1-y0)) - (1-y0);
							p2 = (-1-y0)*log(sqrt((1-x0)^2 + (-1-y0)^2) + (1-x0)) + (1-x0)*log(sqrt((1-x0)^2 + (-1-y0)^2) + (-1-y0)) - (-1-y0);
							p3 = (1-y0)*log(sqrt((-1-x0)^2 + (1-y0)^2) + (-1-x0)) + (-1-x0)* log(sqrt((-1-x0)^2 + (1-y0)^2) + (1-y0)) - (1-y0);
							p4 = (-1-y0)*log(sqrt((-1-x0)^2+(-1-y0)^2)+(-1-x0)) + (-1-x0)* log(sqrt((-1-x0)^2 + (-1-y0)^2) + (-1-y0)) - (-1-y0);
						else
							p1 = -(z-z0)*atan(((1-x0)*(1-y0))/((z-z0)*sqrt(((z-z0)^2 + (1-x0)^2) + (1-y0)^2))) + (z-z0)*atan((1-y0)/(z-z0)) + (1-y0)*log((1-x0) + sqrt((z-z0)^2 + (1-x0)^2 + (1-y0)^2)) + (1-x0)*log(sqrt((z-z0)^2 + (1-x0)^2 + (1-y0)^2) + (1-y0)) - (1-y0);
							p2 = -(z-z0)*atan(((1-x0)*(-1-y0))/((z-z0)*sqrt(((z-z0)^2 + (1-x0)^2) + (-1-y0)^2))) + (z-z0)*atan((-1-y0)/(z-z0)) + (-1-y0)*log((1-x0) + sqrt((z-z0)^2 + (1-x0)^2 + (-1-y0)^2)) + (1-x0)*log(sqrt((z-z0)^2 + (1-x0)^2 + (-1-y0)^2) + (-1-y0)) - (-1-y0);
							p3 =	-(z-z0) *atan(((-1-x0)*(1-y0))/((z-z0)*sqrt(((z-z0)^2 + (-1-x0)^2) + (1-y0)^2))) + (z-z0)*atan((1-y0)/(z-z0)) + (1-y0)*log((-1-x0) + sqrt((z-z0)^2 + (-1-x0)^2 + (1-y0)^2)) + (-1-x0)*log(sqrt((z-z0)^2 + (-1-x0)^2 + (1-y0)^2) + (1-y0)) - (1-y0);
							p4 = -(z-z0) *atan(((-1-x0)*(-1-y0))/((z-z0)*sqrt(((z-z0)^2 + (-1-x0)^2) + (-1-y0)^2))) + (z-z0)*atan((-1-y0)/(z-z0)) + (-1-y0)*log((-1-x0) + sqrt((z-z0)^2 + (-1-x0)^2 + (-1-y0)^2)) + (-1-x0)*log(sqrt((z-z0)^2 + (-1-x0)^2 + (-1-y0)^2) + (-1-y0)) - (-1-y0);
						end;
					else
						if ((z~=z0))
							p1 = -(z-z0)*atan(((1-x0)*(1-y0))/((z-z0)*sqrt(((z-z0)^2 + (1-x0)^2) + (1-y0)^2))) + (z-z0)*atan((1-y0)/(z-z0)) + (1-y0)*log((1-x0) + sqrt((z-z0)^2 + (1-x0)^2 + (1-y0)^2)) + (1-x0)*log(sqrt((z-z0)^2 + (1-x0)^2 + (1-y0)^2) + (1-y0)) - (1-y0);
							p2 = -(z-z0)*atan(((1-x0)*(-1-y0))/((z-z0)*sqrt(((z-z0)^2 + (1-x0)^2) + (-1-y0)^2))) + (z-z0)*atan((-1-y0)/(z-z0)) + (-1-y0)*log((1-x0) + sqrt((z-z0)^2 + (1-x0)^2 + (-1-y0)^2)) + (1-x0)*log(sqrt((z-z0)^2 + (1-x0)^2 + (-1-y0)^2) + (-1-y0)) - (-1-y0);
							p3 =	-(z-z0) *atan(((-1-x0)*(1-y0))/((z-z0)*sqrt(((z-z0)^2 + (-1-x0)^2) + (1-y0)^2))) + (z-z0)*atan((1-y0)/(z-z0)) + (1-y0)*log((-1-x0) + sqrt((z-z0)^2 + (-1-x0)^2 + (1-y0)^2)) + (-1-x0)*log(sqrt((z-z0)^2 + (-1-x0)^2 + (1-y0)^2) + (1-y0)) - (1-y0);
							p4 = -(z-z0) *atan(((-1-x0)*(-1-y0))/((z-z0)*sqrt(((z-z0)^2 + (-1-x0)^2) + (-1-y0)^2))) + (z-z0)*atan((-1-y0)/(z-z0)) + (-1-y0)*log((-1-x0) + sqrt((z-z0)^2 + (-1-x0)^2 + (-1-y0)^2)) + (-1-x0)*log(sqrt((z-z0)^2 + (-1-x0)^2 + (-1-y0)^2) + (-1-y0)) - (-1-y0);
						else % z == z0
							if ( (abs(x0) == 1) && (abs(y0) ~= 1))
						%		log(sqrt((y-y0)^2)) - log(sqrt((y-y0)^2 + 2^2) + 2)
								p1 = -((1/2)*(1-y0)*(log((1-y0)^2) - 2));
								p2 = -((1-y0)*(log(sqrt((1-y0)^2 + 4) + 2) - 1) + 2*asinh((1-y0)/2));
								p3 = -((1/2) *(-1-y0)*(log((-1-y0)^2) - 2));
								p4 = -((-1-y0)*(log(sqrt((-1-y0)^2 + 4) + 2) - 1) + 2*asinh((-1-y0)/2));
							end;
							if ( (abs(x0) ~= 1) && (abs(y0) == 1))
								r0 = x0; % swap them
								x0 = y0;
								y0 = r0;
								p1 = -((1/2)*(1-y0)*(log((1-y0)^2) - 2));
								p2 = -((1-y0)*(log(sqrt((1-y0)^2 + 4) + 2) - 1) + 2*asinh((1-y0)/2));
								p3 = -((1/2) *(-1-y0)*(log((-1-y0)^2) - 2));
								p4 = -((-1-y0)*(log(sqrt((-1-y0)^2 + 4) + 2) - 1) + 2*asinh((-1-y0)/2));
								r0 = x0; % unswap them
								x0 = y0;
								y0 = r0;
							end;
						end;
					end;
				end;
				if (flag == 1)
					p1=0;
					p2=0;
					p3=0;
					p4=0;
					constij(i,j) = 4*asinh(1);
				else
					constij(i,j) = p1 - p2 - p3 + p4;
				end;
			else
				constij(i,j) = 0.0;
			end;
			% TESTED

			% then do the xx term % integrate xx/(x-x0^^2 + (y-y0)^2 + (z-z0)^2)^(3/2) dx dy
			if (i==j) && (i ==ndir)
			% integrate from x = -1 to x = 1 integrate from y = -1 to y = 1  (1/(x^2+y^2+z^2)^(3/2))dy dx
			% integrate 1/(x^2+R^2)^(3/2)dx =  x/(R^2 sqrt(R^2 + x^2))
			% here R^2 = ((y-y0)^2 + (z-z0)^2) and x = (1-x0) or (-1-x0)
      % (1-x0)/( ((y-y0)^2 + (z-z0)^2)*sqrt((y-y0)^2 + (z-z0)^2 + (1-x0)^2)) +
			% (1+x0)/( ((y-y0)^2 + (z-z0)^2)*sqrt((y-y0)^2 + (z-z0)^2 + (1+x0)^2))
			% now integrate in y from -1 to 1
			% integrate 1/( (y^2 + S^2)*sqrt(y^2 + R^2)), if R ~= S
			%  (atan((y sqrt(R^2 - S^2))/(S sqrt(R^2 + y^2))))/(S sqrt(R^2 - S^2)
			% here y = (1-y0) or (-1-y0), S =(z-z0) and R^2 = ((z-z0)^2 + (1-x0)^2)
				
				indx0 = mod(ndir+1,3);
				if (indx0 == 0)
					indx0 = 3;
				end;
				indy0 = mod(ndir+2,3);
				if (indy0 == 0)
					indy0 = 3;
				end;
				x0 = pos(indx0); % not ndir
				y0 = pos(indy0); % not ndir, or x0
				z0 = pos(ndir);
				z = 0; % because we took the difference already between the 2 vectors
				
				% if z==z0, it does not matter, since you will multiply by 0 anyway
									
				if (z~=z0)
					if (x0 ~= 1)
						px1 = (1-x0)*(atan(((1-y0)*(1-x0))/((z-z0)*sqrt((z-z0)^2 + (1-x0)^2 + (1-y0)^2))))/((z-z0)*(1-x0));
						px2 = (1-x0)*(atan(((1+y0)*(1-x0))/((z-z0)*sqrt((z-z0)^2 + (1-x0)^2 + (1+y0)^2))))/((z-z0)*(1-x0));
					else
						px1 = 0; % (1-x0)*(1-y0)/((z-z0)^2*sqrt((z-z0)^2 + (1-y0)^2));
						px2 = 0; % (1-x0)*(1+y0)/((z-z0)^2*sqrt((z-z0)^2 + (1+y0)^2));
					end;
					if (x0 ~= -1)
						px3 = (1+x0)*(atan(((1-y0)*(1+x0))/((z-z0)*sqrt((z-z0)^2 + (1+x0)^2 + (1-y0)^2))))/((z-z0)*(1+x0));
						px4 = (1+x0)*(atan(((1+y0)*(1+x0))/((z-z0)*sqrt((z-z0)^2 + (1+x0)^2 + (1+y0)^2))))/((z-z0)*(1+x0));
					else
						px3 = 0; % (1+x0)*(1-y0)/((z-z0)^2*sqrt((z-z0)^2 + (1-y0)^2));
						px4 = 0; % (1+x0)*(1+y0)/((z-z0)^2*sqrt((z-z0)^2 + (1+y0)^2));
					end;
				else % does not matter since this will be multiplied by zero
					px1 = 0;
					px2 = 0;
					px3 = 0;
					px4 = 0;
				end;
																										
				xxij(i,j) = (z-z0)^2*(px1 + px2 + px3 +px4);
			
% TESTED
			end;
			
			if (i==j) && (i~=ndir)
			% integrate from x = -1 to x = 1 integrate from y = -1 to y = x^2 (1/(x^2+y^2)^(3/2))dy dx
			% SET xo, y0, and z0 properly
			% integrate x^2/(x^2+R^2)^(3/2)dx =   log(sqrt(R^2 + x^2) + x) - x/sqrt(R^2 + x^2)
			% here R^2 = ((y-y0)^2 + (z-z0)^2) and x = (1-x0) or (-1-x0) (1+x0)
      % log(sqrt((y-y0)^2 + (z-z0)^2 + (1-x0)^2) + (1-x0)) - (1-x0)/sqrt((y-y0)^2 + (z-z0)^2 + (1-x0)^2)
			% - log(sqrt((y-y0)^2 + (z-z0)^2 + (1+x0)^2) + (-1-x0)) - (-1-x0)/sqrt((y-y0)^2 + (z-z0)^2 + (1+x0)^2)
			% now integrate in y from -1 to 1
			% integrate log(sqrt(y^2 + R^2) + S)  - S/sqrt(y^2 + R^2) dy =
			% -sqrt(R^2 - S^2) tan^(-1)((S y)/(sqrt(R^2 - S^2) sqrt(R^2 + y^2))) + sqrt(R^2 - S^2) tan^(-1)(y/sqrt(R^2 - S^2)) + y (log(sqrt(R^2 + y^2) + S) - 1)
			% here R^2 = (z-z0)^2 + (1-x0)^2, S = (1-x0) and y = (1-y0) or (-1-y0)
			% and if R== S so if z == z0, but I do not need it unless abs(x0) == 1, apparently, when R == 0 then we get
			%  1/2 y (log(y^2) - 2)
				x0 = pos(i); % the one we have in the numerator
				y0 = pos(6-ndir-i); % the other one
				z0 = pos(ndir);
				z = 0; % because we took the difference already between the 2 vectors
				
				if (z~=z0)
					px1 = -(z-z0)*atan(((1-x0)*(1-y0))/((z-z0)*sqrt((z-z0)^2 +(1-x0)^2 + (1-y0)^2))) + (z-z0)*atan((1-y0)/(z-z0)) + (1-y0)*(log(sqrt((z-z0)^2 + (1-x0)^2 + (1-y0)^2) +  (1-x0)) - 1);
					px2 = -(z-z0)*atan(((1-x0)*(-1-y0))/((z-z0)*sqrt((z-z0)^2 + (1-x0)^2 + (-1-y0)^2))) + (z-z0)*atan((-1-y0)/(z-z0)) + (-1-y0)*(log(sqrt((z-z0)^2 + (1-x0)^2 + (-1-y0)^2) +  (1-x0)) - 1);
					px3 = -(z-z0)*atan(((-1-x0)*(1-y0))/((z-z0)*sqrt((z-z0)^2 + (-1-x0)^2 + (1-y0)^2))) + (z-z0)*atan((1-y0)/(z-z0)) + (1-y0)*(log(sqrt((z-z0)^2 + (-1-x0)^2 + (1-y0)^2) +  (-1-x0)) - 1);
					px4 = -(z-z0)*atan(((-1-x0)*(-1-y0))/((z-z0)*sqrt((z-z0)^2 + (-1-x0)^2 + (-1-y0)^2))) + (z-z0)*atan((-1-y0)/(z-z0)) + (-1-y0)*(log(sqrt((z-z0)^2 + (-1-x0)^2 + (-1-y0)^2) +  (-1-x0)) - 1);
				else % here z == z0
					if (x0 == 1)
						if (y0 ~= 1)
							px1 = 1/2*(1-y0)*(log((1-y0)^2) - 2);
						else
							px1 = 0;
						end;
						if (y0 ~= -1)
							px2 = 1/2*(-1-y0)*(log((-1-y0)^2) - 2);
						else
							px2 = 0;
						end;
					else
						if (y0 ~= 1)
							px1 = (1-y0)*(log(sqrt( (1-x0)^2 + (1-y0)^2) +  (1-x0)) - 1);
						else
							px1 = 0;
						end;
						if (y0 ~= -1)
							px2 =  (-1-y0)*(log(sqrt((1-x0)^2 + (-1-y0)^2) +  (1-x0)) - 1);
						else
							px2 = 0;
						end;
					end;
					if (x0 == -1)
						if (y0 ~= 1)
							px3 = (1-y0)*(log(sqrt((-1-x0)^2 + (1-y0)^2) +  (-1-x0)) - 1);
						else
							px3 = 0;
						end;
						if (y0 ~= -1)
							px4 = (-1-y0)*(log(sqrt((-1-x0)^2 + (-1-y0)^2) +  (-1-x0)) - 1);
						else
							px4 = 0;
						end;
					else
						if (y0 ~= 1)
							px3 = (1-y0)*(log(sqrt((-1-x0)^2 + (1-y0)^2) +  (-1-x0)) - 1);
						else
							px3 = 0;
						end;
						if (y0 ~= -1)
							px4 = (-1-y0)*(log(sqrt((-1-x0)^2 + (-1-y0)^2) +  (-1-x0)) - 1);
						else
							px4 = 0;
						end;
					end;
				end;
				% only problem is if y0 == 1 or y0 == -1
				% TESTED
				xxij(i,j) = (px1 - px2 - px3 +px4);

			end;
			
			if ((i~=j) && (i~=ndir) && (j == ndir))
			% integrate from x = -1 to x = 1 integrate from y = -1 to y = 1 (x/(x^2+y^2)^(3/2))dy dx
			% integrate (x/(x^2+R^2)^(3/2)) dx =  -1/sqrt(R^2 + x^2)
			% here R^2 = (y-y0)^2 + (z-z0)^2 and x is (x-x0)
			% -1/sqrt((y-y0)^2 + (z-z0)^2 + (x-x0)^2)
			% -1/sqrt((y-y0)^2 + (z-z0)^2 + (1-x0)^2)  + 1/sqrt((y-y0)^2 + (z-z0)^2 + (1+x0)^2)
			% integrate (-1/sqrt(y^2 + R^2)) dy =  -log(sqrt(R^2 + y^2) + y)
			% here R^2 = (z-z0)^2 + (1-x0)^2  or R^2 = (z-z0)^2 + (1+x0)^2
				
				x0 = pos(i); % the one we have in the numerator
				y0 = pos(6-ndir-i); % the other one
				z0 = pos(ndir);
				z = 0; % because we took the difference already between the 2 vectors
				
				if (z~=z0)
					px1 =  -log(sqrt((z-z0)^2 + (1-x0)^2  + (1-y0)^2) + (1-y0));
					px2 =  log(sqrt((z-z0)^2 + (1+x0)^2  + (1-y0)^2) + (1-y0));
					px3 =  log(sqrt((z-z0)^2 + (1-x0)^2  + (-1-y0)^2) + (-1-y0));
					px4 =  -log(sqrt((z-z0)^2 + (1+x0)^2  + (-1-y0)^2) + (-1-y0));
				else
					px1 = 0; % not really, but it gets multiplied by 0
					px2 = 0;
					px3 = 0;
					px4 = 0;
				end;

				xxij(i,j) = (z-z0)*(px1 + px2 + px3 +px4);
				% TESTED
			
			end;
			if ((i~=j) && (j~=ndir) && (i == ndir)) % integrate from x = -1 to x = 1 integrate from y = -1 to y = 1 (x/(x^2+y^2)^(3/2))dy dx
				x0 = pos(j); % the one we have in the numerator
				y0 = pos(6-ndir-j); % the other one
				z0 = pos(ndir);
				z = 0; % because we took the difference already between the 2 vectors
				
				if (z~=z0)
					px1 =  -log(sqrt((z-z0)^2 + (1-x0)^2  + (1-y0)^2) + (1-y0));
					px2 =  log(sqrt((z-z0)^2 + (1+x0)^2  + (1-y0)^2) + (1-y0));
					px3 =  log(sqrt((z-z0)^2 + (1-x0)^2  + (-1-y0)^2) + (-1-y0));
					px4 =  -log(sqrt((z-z0)^2 + (1+x0)^2  + (-1-y0)^2) + (-1-y0));
				else
					px1 = 0; % not really, but it gets multiplied by 0
					px2 = 0;
					px3 = 0;
					px4 = 0;
				end;
				xxij(i,j) = (z-z0)*(px1 + px2 + px3 +px4);
				% TESTED
			end;
			if ((i~=j) && (j~=ndir) && (i~=ndir)) % integrate from x = -1 to x = 1 integrate from y = -1 to y = 1  (xy/(x^2+y^2+z^2)^(3/2))dy dx
      % integrate (yx/(x^2+R^2)^(3/2)) dx =  -y/sqrt(R^2 + x^2)
			% here R^2 = (y-y0)^2 + (z-z0)^2 and x is (x-x0), y is (y-y0)
			% -(y-y0)/sqrt((y-y0)^2 + (z-z0)^2 + (x-x0)^2)
			% -(y-y0)/sqrt((y-y0)^2 + (z-z0)^2 + (1-x0)^2)  + (y-y0)/sqrt((y-y0)^2 + (z-z0)^2 + (1+x0)^2)
			% integrate (-y/sqrt(y^2 + R^2)) dy =  -sqrt(R^2 + y^2)
			% here R^2 = (z-z0)^2 + (1-x0)^2  or R^2 = (z-z0)^2 + (1+x0)^2
				x0 = pos(i); % the one we have in the numerator
				y0 = pos(j); % the other one
				z0 = pos(ndir);
				z = 0; % because we took the difference already between the 2 vectors
				
			  px1 = -sqrt((z-z0)^2 + (1-x0)^2 + (1-y0)^2);
				px2 = sqrt((z-z0)^2 + (1-x0)^2 + (1+y0)^2);
				px3 = sqrt((z-z0)^2 + (1+x0)^2 + (1-y0)^2);
				px4 = -sqrt((z-z0)^2 + (1+x0)^2 + (1+y0)^2);
				xxij(i,j) = (px1 + px2 + px3 +px4);
				% TESTED
			end;
			
		end;

	end;
end;

