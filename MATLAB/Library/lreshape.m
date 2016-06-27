function [T,dl,dr]=lreshape(T,llegs,rlegs)
% [T,dl,dr]=lreshape(T,llegs,rlegs)
% Reshape the llegs into one leg and the rlegs into one leg
% Output: reshaped tensor and dimensions of legs
% perform a permutation if needed

%normal tensor
nl = numel(llegs);
nr = numel(rlegs);
nlegs = nl + nr;
%permute legs first
T = permute(T,[llegs rlegs]);
s = size(T);
%adapt size in case of dim 1 legs
s = [s ones(1,nlegs-numel(s))];
dl = s(1:nl);
dr = s(nl+1:nlegs); 
T=reshape(T,prod(dl),prod(dr));

end