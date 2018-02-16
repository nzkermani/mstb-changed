function [ fInfo,numV,tStamp ] = imzmlFileInfo(fold)
%imzmlFileInfo - determine the number of variables and perhaps other
%information in a bunch of processed mat files for the toolbox

% Which folder?
if isempty(fold)
    fold = '/Volumes/JSM/DB/PBC/imzML/';
end

% Get all of the files within
allF = fileFinderAll(fold,'imzML',true);
numF = size(allF,1);

% Save information
fInfo = cell(numF,3);
numV = zeros(numF,1);
tStamp = zeros(numF,1);

% Options
opts.mzFrac   = 0.01;
opts.mzRes    = 0.002;
opts.ppmRes   = 8;
opts.mzRange  = [100 1000];
opts.nazanin  = 0;
opts.prompt = {'Fraction of empty pixels per m/z',...
    'm/z resolution',...
    'Maximum ±ppm shift',...
    'm/z range',...
    'Nazanin? (1 or 0)'};


% Loop through...
for n = 1:numF    
    
    % File information
    file.dir = [allF{n,1} filesep];
    file.nam = allF{n,2};
    
    % How many variables
    %[mz,~,tStamp(n,1)] = desiFunctionIMZML(file,opts);
    mz = 1;
    tStamp(n,1) = imzmlDate([file.dir file.nam]);
    numV(n,1) = numel(mz);    
    
    % File information
    fInfo{n,1} = allF{n,2};
    [~,fInfo{n,2}] = previousFolder(allF{n,1});
    
    disp(int2str(n));
    
end


end

