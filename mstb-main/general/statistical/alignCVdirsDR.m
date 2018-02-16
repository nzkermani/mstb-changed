function [dirIndcs,dirSigns] = alignCVdirsDR(weights1,weights2)

nComps   = size(weights1,2); 
indcs    = 1:nComps;
dirIndcs = 1:nComps;
dirSigns = ones(1,nComps);

for i = 1:nComps
    cossvals    = weights1(:,i)'*weights2;
    [val,idx]   = max(abs(cossvals));
    dirIndcs(i) = indcs(idx(1));
    if cossvals(idx)<0
        dirSigns(i) = -1;
    end 
    indcs(idx(1))  = [];      
    weights2(:,idx) = [];
end
return;