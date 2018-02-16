function [ output_args ] = resizeImages(fold,type,factor)
%resizeImages - read in images from a folder, resize them, and export to
%separate folder...

% Add filesep
if strcmp(fold(end),filesep)
    fold = [fold filesep];
end
if ~exist(fold)
    disp('Does not exist');
    return
end

% Get the files
allF = fileFinderAll(fold,type,true);
numF = size(allF,1);

% Storage location...
newDir = [fold 'LowRes/'];
if ~exist(newDir)
    mkdir(newDir);
end

% Loop through
for n = 1:numF
    
    im = imread([allF{n,1} filesep allF{n,2}]);
    
    im2 = imresize(im,factor);
    
    imwrite(im2,[newDir allF{n,2}]);
    
    % What is the name of this number?
    thisNo = str2double(allF{n,2}(1:end-3));    
    imwrite(im2,[newDir int2str((2*numF)-thisNo+1) '.jpg']);
    
end
    
    


end

