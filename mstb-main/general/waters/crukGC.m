%% Cancer Research Grand Challenge
% Workflow for processing a bunch of Waters Xevo files.

% The raw Waters files have been converted into H5 files in order to make
% them readable on a Mac computer. These files are then processed en masse
% in order to accurately determine the common peaks and then to perform
% multivariate statistics on the data.

%% Define the folder of H5 files
if ismac
    %fold = '/Users/jmckenzi/Documents/Box Sync/H5 Convert/Neg/';
    fold = '/Users/jmckenzi/Documents/Box Sync/H5 Convert/Ovarian-Neg/';
    %fold = '/Volumes/Data/Data/Admin/DESI-Raman Dump/';
    %fold = '/Users/jmckenzi/Desktop/DESIRAMANH5TEST/';
else
    fold = 'E:\Box Sync\H5 Convert\Neg\';
    %fold = 'E:\Data\DESI Raman\H5 Convert\Pos\';
end

% Over what m/z range do we want to analyse?
mzRange = [600 1000];

%% First we need to determine the mean spectrum of each file
[spec,times1] = rawH5global(fold,mzRange);

%% Now resample the m/z axes so that they are all the same
[resamp] = rawH5globalPick(spec);

%% Now we need to align the peaks and perform peak picking
[picked] = rawH5align(resamp);

%% Pick peaks from each file
rawH5extractBatch(fold,spec(:,2),picked);

%% Import from H5
%[gc] = procH5import;

%% Run MVA
%[mva] = procH5mva(gc,'pca-in');

%% Tile results
%tileH5mva(gc,mva);