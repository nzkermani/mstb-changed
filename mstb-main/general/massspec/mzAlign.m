function [data,cmz,qcmz] = mzAlign(data,ppmTol,plotFig)
% First define the cmz vector, and then align everything to it

if nargin == 2
    plotFig = true;
else
    plotFig = false;
end

% Define the peak matching options
opts.ppmTol = ppmTol;
%opts.mzRes = mzRes;
opts.estimationMethod = 'histogram';
opts.display = false;
opts.handBand = eval(['@(z) (' num2str(opts.ppmTol) '* z / 1e6)']);            
opts.handGap  = eval(['@(z) (' num2str(opts.ppmTol) '* z(:,1) / 1e6)']);
opts.handDist = @(R,S) abs(sum((R-S),2));
opts.boolSC = false; 
opts.boolSA = false; 
opts.boolSN = false;

% Prep the mz vectors
numF = size(data,2);
if numF > 1
    mzPeaks = cell(numF,1);
    for n = 1:numF
        try
            mzPeaks{n,1} = [data(n).mz' nanmean(data(n).sp,1)'];
        catch
            mzPeaks{n,1} = [data(n).mz nanmean(data(n).sp,1)'];
        end
    end
    
    % Determine the cmz vector in the traditional way
    [cmz] = mzCMZ(data,opts.ppmTol / 2);

elseif numF == 1
    numF = size(data,1);
    mzPeaks = data;
    
    % Determine the cmz vector
    [cmz] = mzCMZ(mzPeaks,opts.ppmTol / 2);
end

% Add a vector of ones to make it work
cmz = [cmz ones(numel(cmz),1)];

% Vector in which to save the 'original' m/z values of matched peaks.  It
% will be good for QC
qcmz = nan(numF,size(cmz,1));

% Run the sample align function for the first time in order to match the
% peaks and determine the qcmz matrix. I don't particularly care about the
% spectral matrix this time; instead I want to know about the m/z values
% matched to each peak in the cmz.
for n = 1:numF
    
    % Simple alignment
    [j,k] = samplealign2(...
        cmz,...
        mzPeaks{n,1},...
        'Band',opts.handBand,...
        'Gap',opts.handGap,...
        'Distance',opts.handDist,...
        'Quantile',[],...
        'SHOWCONSTRAINTS',opts.boolSC,...
        'SHOWNETWORK',opts.boolSN,...
        'SHOWALIGNMENT',opts.boolSA);

    % qcmz values?
    qcmz(n,j) = data(n).mz(k);

end

% The values from cmz (above) are not accurate, as they are determined via
% an increment in ppm distance from one m/z to the next. Thus they do not
% correspond to any measured value, and may shift somewhat from the raw m/z
% values. Thus here we determine the median m/z value for each
newMZ = nanmean(qcmz,1)';

% Here we run the same function but now we use the better values of mz
for n = 1:numF
    
    [j,k] = samplealign2(...
        [newMZ cmz(:,2)],...
        mzPeaks{n,1},...
        'Band',opts.handBand,...
        'Gap',opts.handGap,...
        'Distance',opts.handDist,...
        'Quantile',[],...
        'SHOWCONSTRAINTS',opts.boolSC,...
        'SHOWNETWORK',opts.boolSN,...
        'SHOWALIGNMENT',opts.boolSA);

    % Make the matrices to finish...
    newSP = zeros(size(data(n).sp,1),size(cmz,1));
    newSP(:,j) = data(n).sp(:,k);
    
    % qcmz values?
    qcmz(n,j) = data(n).mz(k);
    
    % Save to the structure
    data(n).al = newSP;

end

cmz = newMZ;

% Now a subfunction for the alignment figure

if ~plotFig
    return
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alignFig(data)

% Need to determine if we have dates for all of the files... Can't do it if
% only a few do. So gather the date from the files...
numF = size(data,2);
if isfield(data(1),'date')
    dateVals = vertvat(data.date);
    if numel(dateVals) == numF
        yVals = [dateVals; max(dateVals) +7];
        dateFlag = true;
    else
        yVals = 1:numF+1;
        dateFlag = false;
    end
else
    yVals = 1:numF+1;
    dateFlag = false;
end

% QC plot 2
figure; hold on;

for n = 1:numF
    
    % Spectra
    x = mzPeaks{n}(:,1);
    y = mzPeaks{n}(:,2);
    y = y / max(y);
    cidx = round(y * 9) + 1;

    % THis is the yaxis data
    if dateFlag
        y2 = repmat(data(n).date,[numel(x) 1]);
        yVals(1,n) = data(n).date;
    else
        y2 = repmat(n,[numel(x) 1]);
    end
        
    
    scatter(x,y2,120,cidx,'o');
    
end
% Add a couple of days and make the last date value
if dateFlag
    yVals(1,numF+1) = max(yVals) + 7;
end

% Scatter plot the cmz
if dateFlag
    scatter(cmz,repmat(yVals(end),[numel(cmz) 1]),150,'k','d','filled');
else
    scatter(cmz,repmat(n+1,[numel(cmz) 1]),150,'k','d','filled')
end


%yval = 1:numF+1;
for n = 1:size(qcmz,2)
    
    fx = qcmz(:,n) > 0;
    tmp = qcmz(fx,n);
    
    if numel(tmp) > 0
        
        xx = [tmp; cmz(n)];
        yy = [yVals(fx) yVals(end)];
        
        [~,ss] = sort(yy);
        
        plot(xx(ss),yy(ss),'-or',...
            'MarkerFaceColor','red');
        
    end

end

if dateFlag
    datetick('y','yy-mm');
end




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%