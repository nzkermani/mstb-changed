function [mzNew,SpNew] = alignSpToCommonMZ(mz,Sp,mzres)

% align peaks to common range
mzNew = floor(mz.*(1/mzres))./(1/mzres);
mz    = floor(mz.*(1/mzres))./(1/mzres);
mzNew = sort(unique(mzNew))';

mzNew(any(mzNew,2)) = [];
[nSmpls,~]          = size(Sp);
SpNew               = zeros(nSmpls,length(mzNew));

for i = 1:nSmpls
    [~,indcsmzNew,indcsmz] = intersect(mzNew,mz(i,:));
    SpNew(i,indcsmzNew)    = Sp(i,indcsmz);
end

end