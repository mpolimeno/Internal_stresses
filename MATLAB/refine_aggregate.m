function [xc,NC] = refine_aggregate(xcin, NCin);
% This function takes a give cluster and creates a cluster with the same shape but returns an aggregate
% where every cube becomes 8 cubes
% this can be applied to any existing aggregate

j=0; % index of refined cluster
for i=1:NCin
	j=j+1;
	xc(j,:) = 2*xcin(i,:);
	j=j+1;
	xc(j,:) = 2*xcin(i,:)+[2 0 0];
	j=j+1;
	xc(j,:) = 2*xcin(i,:)+[0 2 0];
	j=j+1;
	xc(j,:) = 2*xcin(i,:)+[0 0 2];
	j=j+1;
	xc(j,:) = 2*xcin(i,:)+[2 2 0];
	j=j+1;
	xc(j,:) = 2*xcin(i,:)+[2 0 2];
	j=j+1;
	xc(j,:) = 2*xcin(i,:)+[0 2 2];
	j=j+1;
	xc(j,:) = 2*xcin(i,:)+[2 2 2];
end;

NC = j;
