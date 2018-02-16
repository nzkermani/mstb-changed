function [images] = imzmlIonExtractBatch(fold,mz,ppm)
%imzmlIonExtractBatch - extract a series of ions from a folder's worth of
%imzML files. Save the results in a format such that they can be visualised
%with tileH5mva developed as part of the CRUK GC

% Define the folder containing the imzML files
if isempty(fold)
    %fold = '/Volumes/JSM/DB/DESI colorectal/-ve mode/';
    fold = '/Volumes/JSM/DB/Colorectal DESI/JLA imzML Recalibrated/';
end

% Might be worth adding the Java class path for imzML converter here,
% rather than invoking it each time in the following function.

% Get the files from folder
allF = fileFinderAll(fold,'imzML',true);
%allF = allF(1:20,:)
numF = size(allF,1);


% for n  = 1:numF
%     copyfile([allF{n,1} filesep allF{n,2}],['/Volumes/JSM/DB/DESI Colorectal JLA/' allF{n,2}]);
%     disp(int2str(n));
% end
% images = [];
% % 
% return

% How are we going to store the results?
images = cell(numF,5);

% Loop through...
for n = 1:numF
    
    % What is the filename?
    file = [allF{n,1} filesep allF{n,2}];
    disp([int2str(n) '/' int2str(numF) char(10) allF{n,2}]);
    
    % Run the extraction function
    [images{n,1},images{n,2}] = imzmlIonExtract(file,sort(mz),ppm);    
    images{n,3} = sort(mz);
    images{n,4} = file;
    images{n,5} = imzmlDate([allF{n,1} filesep allF{n,2}]);
    
end

end

