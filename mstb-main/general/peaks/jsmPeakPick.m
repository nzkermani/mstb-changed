function [ mz2,sp2,sp3 ] = jsmPeakPick(mz,sp,method)
%jsmPeakPick - perform various peak picking methods (eventually)
%
%
% INPUTs
% mz
% sp
% method - simple so far

% Save the original data
spOrig = sp;

snr = 10;

% Determine reference spectrum
sp = 1000000 * bsxfun(@rdivide,sp,nansum(sp,2));
refSp = nanmean(sp,1);
unsSp = refSp;

% This should probably be smoothed a little
refSp = masmooth(mz,refSp,10,'mean');

% Baseline smoothing?
[base] = movingWindow(refSp,20,'min');
base = movingWindow(base,20,'mean');

switch method
    case 'locmax'
        [mask] = localMaxima(refSp,base,snr);
        
    otherwise
        % there is no otherwise just yet
end

figure('Units','normalized','Position',[0.25 0.25 0.4 0.4]);
hold on;
plot(mz,refSp,'r');
%plot(mz,unsSp,'m');
plot(mz,base,'k');
plot(mz,base*snr,'k');
scatter(mz(mask),refSp(mask),40,'b','o','filled');

box on;
set(gca,'FontSize',16);
xlabel('m/z','FontSize',18,'FontWeight','bold');
ylabel('Intensity (smoothed)','FontSize',18,'FontWeight','bold');

% Now with the logical matrix of peaks to be included, reduce the dataset
mz2 = mz(mask);
sp2 = spOrig(:,mask);
[sp3] = getIntensityMatrix(spOrig,mask);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [locMax] = localMaxima(refSp,base,snr)

% Create left and right slides
sl = [refSp(2:end) Inf];
sr = [Inf refSp(1:end-1)];

locMax = refSp > sl & refSp > sr & refSp > base*snr;%mean(refSp);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new] = getIntensityMatrix(sp,mask)
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