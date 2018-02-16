function [X,binInd] = dpnTOBG(X,nbins,binInd)
% getObjectPixels identifies pixels that differentiates tissue object from background 
% Refernce: N. Otsu, "A threshold selection method from gray-level histogram, 
%           IEEE Trans on System Man Cybernetics 9 (1979), no. 1, 62-66.
% X - grayscale image 
% nbins - number of bins for histogram estimation

if isempty(nbins)
    nbins = 20;
end

[h,hvals] = hist(X(:),nbins);

% calculation of the threshold as described by N. Otsu, A threshold 
% selection method from gray-level histogram, 
% IEEE Trans on System Man Cybernetics 9 (1979), no. 1, 62-66.
L       = length(h);
i       = 1:L;
A       = cumsum(h);
B       = cumsum(h.*i);
u       = B ./ A;
tmp     = (A(L) - A);
v       = (B(L) - B) ./ (tmp + (tmp==0));
F       = A .* (A(L) - A) .* (u-v).^2;

% This is the threshold for the determination...
if isempty(binInd)
    [~,binInd] = max(F);
end

X(X <=hvals(binInd) - (hvals(2)-hvals(1))/2) = 0;
X(X > hvals(binInd) - (hvals(2)-hvals(1))/2) = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
