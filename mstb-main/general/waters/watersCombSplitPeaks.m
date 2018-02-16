function [X,MZ] = watersCombSplitPeaks(X,MZ,mztol)
% Combine split peaks after rounding off mz values
%
% Review with other versions of combSplitPeak1/2

diffmz = diff(MZ);
if nargin < 3
    mztol = min(diffmz); 
end
mztol = mztol + mztol./10;

[mzuniq,mzindcs]     = find(diffmz>mztol);
dubmzindcs           = find(diff(mzindcs)>1);
if ~isempty(dubmzindcs)
    [nrows,ncols,nvrbls] = size(X);
    X                    = reshape(X,nrows*ncols,nvrbls);
    X_temp = X(:,mzindcs);
    for i = dubmzindcs
        X_temp(:,i+1) = sum(X(:,mzindcs(i)+1:mzindcs(i+1)),2);
    end
    X  = reshape(X_temp,nrows,ncols,length(mzindcs));
    MZ = MZ(mzindcs);
end
return;