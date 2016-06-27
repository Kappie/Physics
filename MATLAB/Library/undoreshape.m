function T=undoreshape(T,dl,dr)
% T=undoreshape(T,dl,dr)
% reshape a matrix back to tensor
% dl: dimensions of left legs (output from lreshape)
% dr: dimensions of right legs (output form lreshape)
    T=reshape(T,[dl dr]);
end
