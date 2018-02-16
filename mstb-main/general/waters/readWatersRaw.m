function [MZ,X,xy,xy2D] = readWatersRaw(filename,mztol,mzfrac,wb,peakDetect)
% reads raw mass spectrometry datasets from Waters Raw files
%   Input: filename - the file location e.g. 'E:\Waters centroid data\DESI_Clinical_combinedAFAMM.raw'
%    Input:  .raw  - the MS image data object of a specimen 
%            mztol  - the mass tolerance of measurements, ...
%                      used to merge peaks than are closer than the mass tolerance of measurements
%            mzfrac - the fraction of pixels with a given MS value (0.001) 
%    Output: MZ - mz values 
%            X  - output matrix [rows x columns x mz]
% Kirill A. Veselkov, Imperial College London, 2015
% Keith Richardson, Waters Corporation, 2015
% James McKenzie, Imperial College London, 2016


if nargin < 2
    mztol = 0.001;
end

if nargin < 3
    mzfrac = 0.01; 
end

if nargin < 5
    peakDetect = 1;
end

% Check if the library is already loaded into Matlab environment
libname = 'MassLynxRaw';

try
    if ~libisloaded(libname)
        loadlibrary(libname,'CMassLynxRawReader.h','addheader','MassLynxRawDefs');
    end
catch err
    error('cannot upload waters .raw data import library: currently works under 64 bit windows');
end

p1       = calllib(libname,'newCMassLynxRawReader',filename);  % pointer to the file 
p2       = calllib('MassLynxRaw','newCMassLynxRawScanReader',p1);
nScans   = calllib('MassLynxRaw','getScansInFunction',p2,1);   % number of scans
sp       = cell(1,nScans); xy = single(zeros(2,nScans));
spcounts = zeros(1,nScans);
irobstd  = zeros(1,nScans);  imedian = zeros(1,nScans);

% Upload data
iterNum = 0;
for iScan = 1:nScans    
    nPoints = calllib('MassLynxRaw','getScanSize',p2,1,iScan); % number of data points in a given spectrum
    mz      = single(zeros(1,nPoints));    % pre-allocate array of zeros for mz values
    int     = single(zeros(1,nPoints));    % pre-allocate array of zeros for intensnities
    mzp     = libpointer('singlePtr',mz);  % get a pointer
    intp    = libpointer('singlePtr',int); % get a pointer
    calllib('MassLynxRaw','readSpectrum',p2,1,iScan,mzp,intp); % read a spectrum
    if peakDetect ==1
        peak      = msPeakFinderV2(double(mzp.Value'),double(intp.Value'),'smoothing',true);
        mzp.Value = peak{1}.mzPeaks(:,1);  intp.Value = peak{1}.peakInt(:,1);
        sp{iScan} = [mzp.Value'; intp.Value'];
    else
        sp{iScan} = [mzp.Value; intp.Value];
    end
       
    % delete zero mz or intensity values 
    irobstd(iScan)  = 1.4826*median(abs(intp.Value - median(intp.Value)));
    imedian(iScan)  = median(intp.Value);
    nonzeroIndcs = mzp.Value~=0 & intp.Value>=0;
    sp{iScan}    = sp{iScan}(:,nonzeroIndcs);
    % round off mz values according to the instrument tolerance
    sp{iScan}(1,:) = round(sp{iScan}(1,:)./mztol); 
    [sp{iScan}(1,:),sp{iScan}(2,:)] = remDuplmz( sp{iScan}(1,:), sp{iScan}(2,:));
    xp  = libpointer('singlePtr',0); % get a pointer for x coordinate
    yp  = libpointer('singlePtr',0); % get a pointer for y coordinate
    calllib('MassLynxRaw','getXYCoordinates',p1,1,iScan,xp,yp);% read a spectrum
    xy(:,iScan)  = [xp.Value yp.Value];
    spcounts(iScan) = length(sp{iScan}(1,:)); 
    
    % now free up the memory space 
    if nargin > 4 && mod(iScan,100) == 0
        frac = iScan/(2*nScans);
        waitbar(frac, wb, ...
            ['Reading from Raw Waters File: ' int2str(iScan) '/' int2str(2*nScans)],...
            'FontSize',6);
    end
    clear mzp intp xp yp 
end

% threshold to filter out noisy peaks
peakThr = mean(imedian) + mean(irobstd)*5;

% Filter out m/z variables occuring less than mzfrac in the dataset 
objectmzs  = zeros(1,sum(spcounts));
mzobjindcs = cumsum([1 spcounts]);

for iScan = 1:nScans
    mz = sp{iScan}(1,:); mz(sp{iScan}(2,:)<=peakThr) = 0;
    objectmzs(mzobjindcs(iScan):mzobjindcs(iScan+1)-1) = mz;
    sp{iScan} = sp{iScan}(1:2,mz~=0);
    % normalize profiles
    sp{iScan}(2,:) = sp{iScan}(2,:)./imedian(iScan);
end

objectmzs = objectmzs(objectmzs~=0);
MZ = unique(objectmzs); %% common mz scale
n  = hist(objectmzs,MZ); clear objectmzs;
% save mz ratios, if the peak is present in at least % of the total sample matrix
MZ = MZ(n > nScans*mzfrac); 

xy2D          = get2Dcoord(xy(1,:),xy(2,:));
[nRows,nCols] = size(xy2D);
X             = zeros(nRows,nCols,length(MZ));
iterNum       = nScans;
for x = 1:nRows
    for y = 1:nCols
        if ~isnan(xy2D(x,y))
            mz  = sp{xy2D(x,y)}(1,:);
            % align mz intensities to a common MZ feature vector
            [mzjoint,mzindcs] = intersect(mz,MZ);
            [ignore,MZindcs]  = intersect(MZ,mzjoint);
            X(x,y,MZindcs)    = sp{xy2D(x,y)}(2,mzindcs);
            
            % Update the waitbar again...
            if nargin > 4 && mod(iterNum,100) == 0
                frac = iterNum/(2*nScans);
                waitbar(frac, wb, ...
                    ['Aligning spectra: ' int2str(iterNum) '/' int2str(2*nScans)]);
            end
                    iterNum = iterNum + 1;
        end
    end
end

MZ     = MZ*mztol;
[X,MZ] = combSplittedPeaks(X,MZ,mztol);
medX   = median(X(X~=0));
% scale median value to 40 for consistency
X      = 40*X./medX; 
return