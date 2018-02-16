function [MZ,X,xy,xy2D,numPoints,opts] = h5waters(filename,opts)
% h5waters - convert a Waters RAW file to an H5 file in order to be able to
% batch process a series of the files, e.g. overnight

% Easily check that this isn't a Mac
if ismac
    error('DOES NOT WORK ON MAC');
end

tic

% Define a special case
special = false;

% What about QC images for raffinose - need to have either positive or
% negative ions, extracted.  We can just do both rather than having to
% decide which
raff = [503.1627 527.1583];

% Default options, or user-definable options to be included later
if nargin == 1
    [opts] = getDefaults;
end

% Here we can decide if we want to acquire a region of interest only. This
% is a new option, so need to ensure that the fucntion is compatible with
% the older versions of the toolbox
if ~isfield(opts,'roi')
    opts.roi = false;
    opts.roiList = [];
end

% Now we can do the ROI part...
if opts.roi && isempty(opts.roiList)
    [opts.roiList,~] = watersScans(filename);
end

% Draw the waitbar
wb = waitbar(0,'h5waters - initialising');

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
if opts.roi
    numS = numel(opts.roiList);
    roi = opts.roiList;
else
    numS = calllib('MassLynxRaw','getScansInFunction',p2,1);
    roi = 1:numS;
end

% Variable initialisation
sp = cell(1,numS); 
xy = single(zeros(2,numS));
spcounts = zeros(1,numS);
totPoints = zeros(1,numS);
irobstd  = zeros(1,numS);  
imedian = zeros(1,numS);

% Images to store raffinose images
raffInt = zeros(2,numS);
raffMZ  = zeros(2,numS);

% This is only for specific cases
if special
    rawSp = cell(numS,1);
end

% Read the scans from the RAW file
for n = 1:numS
    
    % Define the scan as i = roi(n). i != n when only when we do the ROI
    % selection
    i = roi(n);
    if i == 0
        sp{n} = [0; 0];
        continue;
    end
    
    % Number of data points in a given spectrum
    nPoints = calllib('MassLynxRaw','getScanSize',p2,1,i);
    totPoints(i,1) = nPoints;
    
    % Read the xy coordinates for this scan
    xp  = libpointer('singlePtr',0);   
    yp  = libpointer('singlePtr',0);    
    calllib('MassLynxRaw','getXYCoordinates',p1,1,i,xp,yp);    
    xy(:,n)  = [xp.Value yp.Value];
    
    % Can't do anything if there is no data in a pixel - this happens to
    % waters data every now and again
    if nPoints == 0
        sp{n} = [0; 0];
        continue;
    end 
        
    % Preallocate variables
    mz  = single(zeros(1,nPoints));    
    int = single(zeros(1,nPoints));
    
    % Read spectrum
    mzp  = libpointer('singlePtr',mz);    
    intp = libpointer('singlePtr',int);    
    calllib('MassLynxRaw','readSpectrum',p2,1,i,mzp,intp);
    
    % Extract raffinose peak - if in centroid mode then find the closest
    % peak. If in profile mode find the largest within a small window +/-
    % 0.3
    switch opts.method
        case 'Centroid'
            try
                df = bsxfun(@minus,mzp.Value',raff);
                [~,idx] = min(abs(df),[],1);
                raffMZ(:,i) = mzp.Value(idx)';
                raffInt(:,i) = intp.Value(idx)';
            catch
                disp('Missing spectrum');
            end
                        
        case 'Profile'
            
    end
    
    % Trim out the variables which lie outside the specified m/z range
    try
        mask = mzp.Value' > min(opts.mzRange) & mzp.Value' < max(opts.mzRange);
        mzp.Value = mzp.Value(mask);
        intp.Value = intp.Value(mask);
    catch
        continue;
    end
    
    % Save the raw spectra
    if special
        rawSp{n,1} = [mzp.Value' intp.Value'];
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
        sp{n} = [mzp.Value'; intp.Value'];
    
    elseif opts.peakDetect == 1 && nPoints < 50
        % Do nothing - is this an option?
        mzp.Value = min(opts.mzRange);
        intp.Value = 0;
        sp{n} = [mzp.Value; intp.Value];
    else        
        % Save the full spectral data
        try
            sp{n} = [mzp.Value; intp.Value];
        catch
            sp{n} = [0; 0];
            continue;
        end
    end
    
    % Determine mean and std of intensity values
    irobstd(n)  = 1.4826*nanmedian(abs(intp.Value - nanmedian(intp.Value)));    
    imedian(n)  = nanmedian(intp.Value);    
    
    % Select only non-zero mz and intensity values
    nonzeroIndcs = mzp.Value > 0 & intp.Value > 0; %(irobstd(i)/2);
    sp{n} = sp{n}(:,nonzeroIndcs);
    
    % Round mz values according to the instrument tolerance
    sp{n}(1,:) = round(sp{n}(1,:) ./ opts.mzRes); 
    
    % Remove duplicate mz values
    [sp{n}(1,:),sp{n}(2,:)] = remDuplmz( sp{n}(1,:), sp{n}(2,:));
    
    % Count the length, i.e. number of variables in this scan
    spcounts(n) = length(sp{n}(1,:));
    
    
    % Update the waitbar
    if mod(n,100) == 0
        frac = n/(numS);
        waitbar(frac, wb, ...
            ['RAW > H5: ' int2str(n) '/' int2str(numS)],...
            'FontSize',10);
    end
    
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
    sp{i}(2,:) = sp{i}(2,:)./imedian(i);
end

% Filter out the values with a zero intensity
objectmzs = objectmzs(objectmzs~=0);

% Determine unique m/z values, against which all pixels can be aligned
MZ = unique(objectmzs);

% See how frequently each m/z value occurs
n  = hist(objectmzs,MZ);
clear objectmzs;

% Remove m/z values that appear in too few of the pixels
% Is this done too early?
MZ = MZ(n > numS*opts.mzFrac); 

% Determine the xy coordinates for the sample
xy2D          = get2Dcoord(xy(1,:),xy(2,:));
[nRows,nCols] = size(xy2D);

% Preallocate the spectral datacube. Can we make this sparse without 
% affecting the rest of the code?
X = zeros(nRows,nCols,length(MZ));
numPoints = zeros(nRows,nCols);

% IMages for storing the QC of raffinose
qcMZ = NaN(nRows,nCols,2);
qcIn = NaN(nRows,nCols,2);

% Loop through each scan and align the m/z vectors
iterNum = numS;
for x = 1:nRows
    for y = 1:nCols
        if ~isnan(xy2D(x,y))
            
            % Extract this scan's m/z vector
            mz  = sp{xy2D(x,y)}(1,:);
            
            % Align mz values to the master MZ list
            [mzjoint,mzindcs] = intersect(mz,MZ);
            [~,MZindcs] = intersect(MZ,mzjoint);
            
            % Place the data in the cube
            X(nRows-x+1,y,MZindcs) = sp{xy2D(x,y)}(2,mzindcs);
            
            % Put the totPoints into an image
            numPoints(nRows-x+1,y) = totPoints(xy2D(x,y),1);
            
            % QC images...
            qcMZ(nRows-x+1,y,1) = raffMZ(1,xy2D(x,y));
            qcMZ(nRows-x+1,y,2) = raffMZ(2,xy2D(x,y));
            
            qcIn(nRows-x+1,y,1) = raffInt(1,xy2D(x,y));
            qcIn(nRows-x+1,y,2) = raffInt(2,xy2D(x,y));
                        
            % Update the waitbar again...
            if nargin > 4 && mod(iterNum,100) == 0
                frac = iterNum/(2*numS);
                waitbar(frac, wb, ...
                    ['Aligning spectra: ' int2str(iterNum) '/' int2str(2*numS)]);
            end
                    iterNum = iterNum + 1;
        end
    end
end

% Save the images in the opts structure
opts.qc.mz = qcMZ;
opts.qc.sp = qcIn;

% Convert the MZ vector to the proper scale
MZ = MZ * opts.mzRes;

%assignin('base','Xlow',X);
%assignin('base','MZlow',MZ);

% Combine peaks that are split, i.e. slight mis-alignment problems. This
% could perhaps benefit from a ppm-based approach
[X,MZ] = watersCombSplitPeaks(X,MZ,opts.mzRes);

% Perhaps here we run the analysis again but using the image
% complementarity function
[MZ,X] = combCompNeighbours(MZ,X,opts.ppmRes);

% Determine the overall median intensity value
%medX = median(X(X~=0));

% Scale median value to 40 for consistency - this is not particularly
% necessary I think
%X = 40 * X ./ medX;

% Now that we have the processed data, we should save it to a file...
%newName = [filename(1:end-3) 'mat'];
%save(newName,'MZ','X','xy','xy2D','filename','opts');
%disp(newName);

% If special then we can just return the raw data to the workspace
if special
    numPoints = rawSp;
end

% Delete the waitbar
delete(wb);

toc

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = getDefaults

opts.mzRes      = 0.01;
opts.mzFrac     = 0.01;
opts.peakDetect = 0;
opts.ppmRes     = 10;
opts.mzRange    = [100 1000];
opts.roi        = false;
opts.roiList    = [];
opts.recal      = true;
opts.method     = 'Centroid';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%