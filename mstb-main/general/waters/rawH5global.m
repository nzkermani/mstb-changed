function [spec,times] = rawH5global(folder,mzRange)
%rawH5global - import the mean spectra for all files and then put on a
%common scale, smooth etc, and then peak pick... This will be used as the
%reference for all subsequent files...

% Define the m/z range - start small...
if nargin == 1
    mzRange = [500 550];
end

% Get all files
files = fileFinderAll(folder,'h5',true);

files = files([1 2 4 5 6],:);

numF = size(files,1);
%numF = 10;

% Cell to store the spectra
spec = cell(numF,2);

% Times
times = zeros(numF,1);

% Loop through each one
for n = 1:numF
    
    tt = tic;
    
    try
        % File name?
        fn = [files{n,1} filesep files{n,2}];
    
        % Import
        [mz,sp] = rawH5mean(fn,mzRange);
    
        % Save
        spec{n,1} = [mz sp];
        spec{n,2} = files{n,2};
    catch err
        err
        disp('FAIL');
        disp(files{n,2});
        spec{n,2} = files{n,2};
    end    
    
    times(n,1) = toc(tt);
end


end

