function [ pks ] = desiRawPeakPick(mz,sp,flt,rp,fig)
%desiRawPeakPick - peak pick, subject to arbitrary threshold

% Determine peak selection threshold
thr = median(flt) * 5;
%thr = median(flt(flt > 0)) * 10;

% Find the peaks in the 
[fx,fy] = findPeaks(flt,thr);

% Now with the list of peaks, need to determine the peak widths...
[hw,pkw,pkwY] = peakWidths(mz,fx,rp);

% Prepare outputs ready to be extracted
pks.mz = mz(fx);
pks.hw = hw;

if nargin == 5
    % Plot the figure
    figure; hold on;
    plot(mz,flt,'k','LineWidth',2);
    plot(mz,sp,'b','LineWidth',2);
    scatter(mz(fy),flt(fy),200,'blue','o');
    scatter(mz(fx),flt(fx),80,'red','o','filled');
    plot(pkw,pkwY,'LineWidth',2,'Color','r');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fx,fy] = findPeaks(flt,thr)
% Find local maxima, and then remove flat-topped peaks

% Slide the vector places to the left and right...
l1 = [flt(2:end) Inf];
l2 = [flt(3:end) Inf Inf];
r1 = [Inf flt(1:end-1)];
r2 = [Inf Inf flt(1:end-2)];

% Find local maxima
fx = flt > l1 & flt > l2 & ...
    flt >= r1 & flt >= r2 & ...
    flt > thr;

fy = (flt == l1 | flt == r1) & ...
    flt > thr;

% Now loop through the fy instances...
lm = find(fy);
%lm2 = lm;
for n = 1:numel(lm)
    
    % Skip ones done already
    if isnan(lm(n))
        continue;
    end
    
    % Differences
    tmp = lm - lm(n);
    
    % Continue until non sequential
    i = n;
    while (tmp(i+1) == tmp(i) + 1) & i < (numel(lm)-1);
        i = i + 1;
    end
    if i == numel(lm)-1
        i = i + 1;
    end    
    fidx = lm(n:i);    
    
    % Determine which of these is to be kept
    keep = round(median(fidx));
    
    % Remove the others...
    fx(setdiff(fidx,keep)) = false;
    fx(keep) = true;
    
    % Ditch
    lm(n:i) = NaN;
    
end



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hw,pws,pwsY] = peakWidths(mz,fx,rp)
% Determine peak widths

% Use the resolution
hw = (mz(fx) / rp) / 2;

% Create a list of peak widths for the plot
pws = [mz(fx)'-hw' mz(fx)'+hw' NaN(sum(fx),1)]';
pws = pws(:);
pwsY = pws;
pwsY(isnan(pws)) = NaN;
pwsY(~isnan(pws)) = -1e6;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
