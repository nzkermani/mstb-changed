function [new] = neighbourIntensityMatrix(sp,mask)
% Look at neighbouring peaks to the mask peaks and take the biggest...

[numO,~] = size(sp);

% Indices of the peaks
mask = find(mask);
numV = numel(mask);

new = zeros(numO,numV);

for n = 1:numO
    
    for r = 1:numV
        
        i = mask(r);
        
        new(n,r) = max(sp(n,i-1:i+1));
        
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%