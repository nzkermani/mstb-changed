function [ output_args ] = copyFileFolder
%copyFileFolder - do stuffto copy a folder, and also a specific file from
%the Z drive to the local drive.

% Where from?
%orig = '/Volumes/Seagate Expansion Drive/PhD DATA/SOTON -ve/';
%orig = '/Volumes/Data/Lab Data/Liver/PBC/H5 files and row data/';
orig = '/Volumes/Data/Lab Data/Endometrial/Olivia/Olivia DESI Raw data/1 scan per sec';

% Where to?
%locn = '/Users/jmckenzi/Desktop/Liam imzML/';
locn = '/Volumes/JSM/DB/PBC/imzML/';

% Which types of files?
type = 'ibd';

% So find all the mat files
allF = fileFinderAll(orig,type,true);

% Find only those beginning with 'S'
fx = ~cellfun(@isempty,cellfun(@min,strfind(allF(:,2),'S'),'UniformOutput',false));
allF = allF(fx,:);
numF = size(allF,1);

% Loop through
for n = 1:numF
    
    % Determine containing folder
    sl = strfind(allF{n,1},filesep);
    cnt = allF{n,1}(sl(end)+1:end);
    
    % New folder...
    newF = [locn cnt filesep];
    if ~exist(newF,'dir')
        mkdir(newF);
    end
    
    % Don't bother if the file already exists
    if exist([newF allF{n,2}],'file')
        continue;
    end
    
    % Now copy the file...
    copyfile([allF{n,1} filesep allF{n,2}],[newF allF{n,2}]);
    
    disp(int2str(n));
    
end


end

