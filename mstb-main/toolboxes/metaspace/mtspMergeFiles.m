function [fName,numA] = mtspMergeFiles
%mtspMergeFiles - combienn the metaspace annotations with the original data
%processing, using the files from the original attempt.  Note that this
%function needs to be changed for each tissue type.

% Location of the (breast) files
origLocn = '/Volumes/JSM/DB/Metaspace/Engine Dump/MTSP-Merge/Neg/Breast/';

% Where are the new annotations?
newLocn = '/Volumes/JSM/DB/Metaspace/NewEngineDump/Convert/';

% Where to save?
saveLocn = '/Volumes/JSM/DB/Metaspace/NewEngineDump/Merge/Breast/';


% Get list of 'old' files
allF = fileFinderAll(origLocn,'mat',true);
numF = size(allF,1);

% New files
newF = fileFinderAll(newLocn,'mat',true);

% Annotation quantity
numA = zeros(numF,2);

% Loop
for n = 1:numF
        
    % Save file name - can skip if done
    tmpSave = [saveLocn newF{n,2}];
    if exist(tmpSave,'file')
        continue;
    end

    % Construct
    tmpF = [allF{n,2}(1:end-4) '-RECAL.mat'];
    disp(tmpF);

    % Find match
    fx = strcmp(newF(:,2),tmpF);
    
    % If just the one file, then open, ditch old, copy, save, etc
    if sum(fx) ~= 1
        disp('Failure to match');
        continue;
    end
    
    % Open the file
    fNew = open([newF{fx,1} filesep newF{fx,2}]);
    
    % Open the original file
    fOrig = open([allF{n,1} filesep allF{n,2}]);
    dpn = fOrig.dpn;
    %clear fOrig
    
    % Replace the bits...
    [dpn2] = mtspMatch(dpn,fNew);
    
    % Save the dpn...
    dpn.mtsp = dpn2.d1;
    
    % Determine the number of annotations
    numA(n,:) = [size(fOrig.dpn.mtsp.sp,3) size(dpn.mtsp.sp,3)];
    
    % Save
    save(tmpSave,'dpn');
    
    
end

fName = allF(:,2);


end

