function [ picked ] = rawH5align(str)
%rawH5align - spectral alignment of the profiles in order to determine a
%consensus m/z vector for peaks

tt = tic;

% We need to provide some reference peaks in order to correctly align the
% data. This can be high intensity peaks from the mean spectrum
mz = str.reMZ;
sp = bsxfun(@rdivide,str.reSp,nansum(str.reSp,2)) * 1000;
ref = nanmedian(sp,1);

% Find the 'best' peaks for alignment
%refPeaks = [601.1312 750.5431 829.3542 885.5489];
refPeaks = [503.1618 539.1384];

% Alignment...
if numel(refPeaks) > 1
    
    % Use Matlab's own function for determining the shifts and for
    % performing alignment
    [spAl,~,shift] = msalignJSM(mz',sp',refPeaks,...
        'Rescaling',false,...
        'ShowPlot',false,...
        'MaxShift',[-0.2 0.2],...
        'WidthOfPulses',0.1);
    
elseif numel(refPeaks) == 1
    % This is simply a shift according to the difference between the
    % expected peak.  Even when using the main raffinose signal at 503 we
    % could avoid this by using other raffinose adducts. But worth
    % implementing for implementation's take
    [spAl,shift] = alignSinglePeak(mz,sp,refPeaks);
    
end

% Determine the median spectrum
ref = nanmedian(spAl,2)';

% Now we need to find the peaks in here
[pks,ref] = rawH5peaks(mz,ref);

% Match the raw m/z values from the data to the selected peaks
[origMZ,origSp] = match2original(mz,sp,pks,shift);

% QC plot to show aligned peaks
qcPlot(mz,sp,ref,origMZ,origSp,pks);

% Another plot...
%diffMZ = bsxfun(@minus,origMZ,pks.mz');
%diffPPM = bsxfun(@rdivide,diffMZ,pks.mz') * 1e6;
%figure; boxplot(diffPPM);

% Prepare the output of things
picked.pks = pks;
picked.origMZ = origMZ;
picked.refSpec = [mz' ref'];
picked.al = spAl;
picked.time = toc(tt);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [al,shift] = alignSinglePeak(mz,sp,ref)
% Perform alignment based on a single peak

% Find the best match for the peak ref in each spectrum.  This will assume
% that the nearest peak is the correct one. Mega shifts cannot be handled
% with this automated one-peak approach
numF = size(sp,1);
peakMZ = zeros(numF,1);
mzD = mz - ref;
for n = 1:numF
    
    % Smooth the spectrum
    sm = sgolayfilt(sp(n,:),7,21);
    
    % Find a local maximum
    lm = sm >= [sm(2:end) Inf] & sm >= [Inf sm(1:end-1)];
    
    % Find the nearest local maximum
    fx = lm .* mzD;
    fx(fx == 0) = Inf;    
    [~,b] = min(abs(fx));    
    peakMZ(n,1) = mz(b);
    
end

% Now interpolate intensities
al = zeros(size(sp));
for n = 1:numF
    al(n,:) = interp1(mz,sp(n,:),mz-(ref-peakMZ(n,1)));
end

% Prepare outputs. We need to transpose al to keep it in line with the
% other function
shift = peakMZ - ref;
al = al';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [origMZ,origSp] = match2original(mz,sp,pks,shift)
% Determine the m/z in the raw data of the peaks, using the simple offset
% provided by the msalignJSM function

% How many files?
numS = numel(shift);

% How many peaks are there?
numP = numel(pks.mz);

% Store the orginal m/z values in here
origMZ = zeros(numP,numS);
origSp = zeros(numP,numS);

% Loop through each file and each peak
for n = 1:numS
    
    for r = 1:numP
        
        % This is the mz of the peak
        peakMZ = pks.mz(r);
        
        % Apply the offset
        localMZ = peakMZ + shift(n);
        
        % Look for the closest peak to localMZ in the average spectrum from
        % this file.  The modified findpeaksJSM function returns both half
        % widths; to prevent using very broad peaks we use the smallest
        % half-width
        span = abs(pks.span(:,r) - pks.mz(r));
        w = min(span);
        
        % Gather the spectrum within this window
        fx = mz > localMZ-w & mz < localMZ+w;
        tmp = [mz(fx)' sp(n,fx)'];
        
        % Instead I want to find the biggest local maximum. SHould I be
        % finding the closest local maximum instead? This is only a problem
        % with poorly resolved low intensity peaks
        fy = tmp(:,2) > [tmp(2:end,2); 0] & tmp(:,2) > [0; tmp(1:end-1,2)];
        fz = double(fy) .* tmp(:,2);
        
        % Find the largest of them here
        [~,idx] = max(fz);
        
        origMZ(r,n) = tmp(idx,1);
        origSp(r,n) = tmp(idx,2);
        
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function qcPlot(mz,sp,ref,origMZ,origSp,pks)

% Create a vector of integers to group and use for sorting
origPk = repmat([1:size(origMZ,1)]',[1 size(origMZ,2)+1])';
origPk = origPk(:);

% Make a plot to show aligment
qcX = [origMZ NaN(size(origMZ,1),1)]';
qcX = qcX(:);

% Make an y vector of intensities in a similar fashion
qcY = [origSp NaN(size(origSp,1),1)]';
qcY = qcY(:);

% Sort according to peak and then intensity
all = sortrows([origPk qcX qcY],[1 3]);

figure; hold on;
cols = parula(size(sp,1));
for n = 1:size(sp,1)
    plot(mz,sp(n,:),...
        'Color',cols(n,:),...
        'LineWidth',1,...
        'LineStyle','-');
end
    
%plot(mz,-ref,'k','LineWidth',4);
plot(all(:,2),all(:,3),'-ko','LineWidth',2);
%scatter(pks.mz,-pks.ints,80,'r','o','filled')


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%