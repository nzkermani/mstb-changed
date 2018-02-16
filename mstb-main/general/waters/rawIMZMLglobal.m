function [mzVec,avg,frq,files,times] = rawIMZMLglobal(folder,mzRange)
%rawIMZMLglobal - do stuff to import and get mean spectra from imzml files,
%and then consider defining the cmz vector and then aligning individual
%scans to that function

% Define the m/z range - start small...
if nargin == 1
    mzRange = [600 900];
end
mzRes = 0.0001;

% Determine an mz vector of ppm increasing values
if mzRes > 1
    mzVec = ppmVector(min(mzRange),max(mzRange),mzRes);
else
    mzVec = (min(mzRange):mzRes:max(mzRange))';
end

% Get all files
files = fileFinderAll(folder,'imzML');
files(1,:) = [];
numF = size(files,1);

% Average spectra
avg = zeros(numel(mzVec),numF);
frq = zeros(numel(mzVec),numF);

% Times
times = zeros(numF,1);

% Loop through each one
for n = 1:numF
    
    tt = tic;
    
    %try
        % File name?
        fn = [files{n,1} filesep files{n,2}];
    
        % Import
        [~,avg(:,n),frq(:,n)] = rawIMZMLmean(fn,mzVec);
    
    %catch err
        %err
        %disp('FAIL');
        %disp(files{n,2});
    %end    
    
    times(n,1) = toc(tt);
end


end

