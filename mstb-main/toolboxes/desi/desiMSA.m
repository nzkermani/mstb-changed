function [data] = desiMSA(data,opts)
%desiMSA - DESI multiple sample analysis
%
% Load files and align their m/z vectors. Then do other stuff like
% normalisation and stuff and then do stuff with the pixels and whatever

% Load the files
if nargin == 0 || isempty(data)
    % Use default files
    %path = '/Volumes/Data/Data/Jocelyn Waters Data/Raman mouse samples mat files/';
    path = '/Users/jmckenzi/Dropbox/Imperial/Projects/DESI Raman Myelin Brain/DESI Data/Updated/';
    [hits] = findFileType(path,'mat');
    
    % Trim the list (negative)
    %hits = hits([1:2:13 16 18 22 24 28 29],:);
    
    % Trim the list (LPC, negative)
    hits = hits(1:5,:);
    
    % Import the files
    [data] = importFiles(hits);
    return

else
    %error('Ask user for files');
end

% Define options for, e.g., peak matching
opts.ppmTol = 50;
opts.mzRes = 0.1;

% Run the alignment
%[data,cmz] = mzAlignment(data,opts);
[cmz,data] = dbBinning(data,1);

% Combine the spectra together, and make vectors so that we know which file
% each spectrum belongs to.
numF = size(data,2);
allInfo = struct('fileID',[],'histID',[]);
for n = 1:numF
    tmp = data(n).file;
    dot = strfind(tmp,'.');
    allInfo(n).fileID = repmat({tmp(1:dot(end)-1)},[size(data(n).al,1) 1]);
    allInfo(n).histID = data(n).histID;
end
full.mz = cmz;
full.sp = vertcat(data.al);
full.meta.histID = vertcat(allInfo.histID);
full.meta.fileID = vertcat(allInfo.fileID);
data = full;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = importFiles(hits)
% Read in the pertinent parts of the files

% Data storage
numF = size(hits,1);
data = struct('path',[],'file',[],'mz',[],'sp',[],'histID',[],'al',[]);

% Loop through the files
for n = 1:numF
    
    % Read in each one
    tmp = open([hits{n,1} filesep hits{n,2}]);
    data(n).path = hits{n,1};
    data(n).file = hits{n,2};
    
    % Determine the histologically annotated pixels
    [isAnno,histID,~] = desiAnnotationExtract(tmp.dpn);
    mask = isAnno ~= 0;
    
    % Save the MS information
    data(n).mz = tmp.dpn.d1.mz;
    data(n).histID = histID(mask);
    
    % Reshape the matrix
    sz = size(tmp.dpn.d1.sp);
    sp = reshape(tmp.dpn.d1.sp,[sz(1)*sz(2) sz(3)]);
    data(n).sp = sp(mask,:);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,cmz] = mzAlignment(data,opts)
% First define the cmz vector, and then align everything to it

% Define the peak matching options
opts.estimationMethod = 'histogram';
opts.display = true;
opts.handBand = eval(['@(z) (' num2str(opts.ppmTol) '* z / 1e6)']);            
opts.handGap  = eval(['@(z) (' num2str(opts.ppmTol) '* z(:,1) / 1e6)']);
opts.handDist = @(R,S) abs(sum((R-S),2));
opts.boolSC = false; 
opts.boolSA = false; 
opts.boolSN = false;

% Prep the mz vectors
numF = size(data,2);
mzPeaks = cell(numF,1);
for n = 1:numF
    mzPeaks{n,1} = [data(n).mz' nanmean(data(n).sp,1)'];
end

% Determine the cmz vector in the traditional way
cmz = mspmatch(mzPeaks,...
    'estimationMethod',opts.estimationMethod,...
    'mzRes',opts.mzRes,...
    'display',opts.display);

% Now run the samplealign2 function to match A to B...
for n = 1:numF
    
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

    % Make the matrices to finish...
    newSP = zeros(size(data(n).sp,1),size(cmz,1));
    newSP(:,j) = data(n).sp(:,k);
    
    % Save to the structure
    data(n).al = newSP;

end

cmz = cmz(:,1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


