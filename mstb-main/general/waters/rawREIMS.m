function [MZ,X,xy,opts] = rawREIMS(filename,opts)
% rawREIMS - process the raw REIMS data files using a similar approach to
% that developed for the DESI data.  Needs to be carefully tested but is a
% good way to get away from using OMB.

% Let's just get the options...
if nargin == 0 && nargout == 1
    [MZ] = getDefaults;
    return
end
    

disp(filename);

% Easily check that this isn't a Mac
if ismac
    error('DOES NOT WORK ON MAC');
end

% Define a special case
special = true;

% Default options, or user-definable options to be included later
if nargin == 1
    [opts] = getDefaults;
end

% Draw the waitbar
wb = waitbar(0,'Raw REIMS - GO GO GO');

% Replace cal function in header?
if opts.recal
    [origRecal] = watersRecal(filename,'');
    if length(origRecal) == 0
        opts.recal = false;
    end
end

% Default initialisation
[p1,p2] = watersPackages(filename);

% Number of scans
numS   = calllib('MassLynxRaw','getScansInFunction',p2,1);

% Variable initialisation
sp = cell(1,numS); 
spcounts = zeros(1,numS);
totPoints = zeros(numS,1);
irobstd  = zeros(1,numS);  
imedian = zeros(1,numS);

% This is only for specific cases
if special
    rawSp = cell(numS,1);
end

% Read the scans from the RAW file
for i = 1:numS
    
    % Number of data points in a given spectrum
    nPoints = calllib('MassLynxRaw','getScanSize',p2,1,i);
    totPoints(i,1) = nPoints;
    
    % Preallocate variables
    mz      = single(zeros(1,nPoints));    
    int     = single(zeros(1,nPoints));
    
    % Read spectrum
    mzp     = libpointer('singlePtr',mz);    
    intp    = libpointer('singlePtr',int);    
    calllib('MassLynxRaw','readSpectrum',p2,1,i,mzp,intp);
    
    % Trim out the variables which lie outside the specified m/z range
    mask = mzp.Value' > min(opts.mzRange) & mzp.Value' < max(opts.mzRange);
    mzp.Value = mzp.Value(mask);
    intp.Value = intp.Value(mask);
    
    % Save the raw spectra
    if special
        rawSp{i,1} = [mzp.Value' intp.Value'];
    end
    
    % Here we do peak detection
    if opts.peakDetect == 1 && nPoints > 50
        
        % This is the peak detection function
        peak = msPeakFinderV2(double(mzp.Value'),...
            double(intp.Value'),...
            'smoothing',true);

        % Save the detected peaks to the matrix
        mzp.Value = peak{1}.mzPeaks(:,1);  
        intp.Value = peak{1}.peakInt(:,1);
        sp{i} = [mzp.Value'; intp.Value'];
    
    elseif opts.peakDetect == 1 && nPoints < 50
        % Do nothing - is this an option?
        mzp.Value = min(opts.mzRange);
        intp.Value = 0;
        sp{i} = [mzp.Value; intp.Value];
    else        
        % Save the full spectral data
        sp{i} = [mzp.Value; intp.Value];
    end
    
    % Determine mean and std of intensity values
    irobstd(i)  = 1.4826*median(abs(intp.Value - median(intp.Value)));    
    imedian(i)  = median(intp.Value);    
    
    % Select only non-zero mz and intensity values
    nonzeroIndcs = mzp.Value > 0 & intp.Value > 0; %(irobstd(i)/2);
    sp{i} = sp{i}(:,nonzeroIndcs);
    
    % Round mz values according to the instrument tolerance
    sp{i}(1,:) = round(sp{i}(1,:) ./ opts.mzRes); 
    
    % Remove duplicate mz values
    [sp{i}(1,:),sp{i}(2,:)] = remDuplmz( sp{i}(1,:), sp{i}(2,:));
    
    % Count the length, i.e. number of variables in this scan
    spcounts(i) = length(sp{i}(1,:));
    
    % Update the waitbar
    frac = i/(numS);
    waitbar(frac, wb, ...
        ['Raw REIMS - ' int2str(i) '/' int2str(numS)],...
        'FontSize',10);
    
    % Clear the memory - though perhaps unnecessary as it gets reused
    clear mzp intp xp yp 
end

% Restore the original cal function in header (if needed)
if opts.recal
    [~] = watersRecal(filename,origRecal);
end

% Threshold to filter out noisy peaks
peakThr = mean(imedian) + mean(irobstd)*5;

% Start to collate all the m/z values throughout the sample
objectmzs  = zeros(1,sum(spcounts));
mzobjindcs = cumsum([1 spcounts]);

for i = 1:numS
    
    % This is the mz vector for this pixel
    mz = sp{i}(1,:);
    
    % Set intensities less than the threshold to 0
    mz(sp{i}(2,:) <= peakThr) = 0;
    
    % Put the m/z values into a single long vector
    objectmzs(mzobjindcs(i):mzobjindcs(i+1)-1) = mz;
    
    % Save only those non-zero values
    sp{i} = sp{i}(1:2,mz~=0);
    
    % Normalize profiles - why do we have to do this?
    sp{i}(2,:) = sp{i}(2,:) ./ imedian(i);
end

% Filter out the values with a zero intensity
objectmzs = objectmzs(objectmzs~=0);

% Determine unique m/z values, against which all pixels can be aligned
MZ = unique(objectmzs);

% Preallocate the spectral matrix. Can we make this sparse without 
% affecting the rest of the code?
X = zeros(numS,length(MZ));

% Loop through each scan and align the m/z vectors
for x = 1:numS
    % Extract this scan's m/z vector
    mz  = sp{x}(1,:);
            
    % Align mz values to the master MZ list
    [mzjoint,mzindcs] = intersect(mz,MZ);
    [~,MZindcs] = intersect(MZ,mzjoint);
            
    % Place the data in the matrix
    X(x,MZindcs) = sp{x}(2,mzindcs);
            
    % Update the waitbar again...
    frac = x/numS;
    waitbar(frac, wb, ...
        ['Aligning spectra: ' int2str(x) '/' int2str(numS)]);
        
end

% Convert the MZ vector to the proper scale
MZ = MZ * opts.mzRes;

% If special then we can just return the raw data to the workspace
if special
    xy = rawSp;
else
    xy = [];
end

% Delete the waitbar
delete(wb);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = getDefaults

opts.mzRes      = 0.01;
%opts.mzFrac     = 0.01;
opts.peakDetect = 1;
%opts.ppmRes     = 10;
opts.mzRange    = [50 1200];
opts.recal      = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%