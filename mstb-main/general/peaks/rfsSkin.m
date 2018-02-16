function [ avgs ] = rfsSkin
%rfsSkin - analysis for Renata's skin samples

% File information
fold = '/Volumes/JSM/DB/Renata/';
file = {'Skin-2C-Full.mat';'Skin-2C-NR5.mat';'Skin-2C-NR15.mat'};
numF = numel(file);

% Save averages
avgs = cell(numF,3);

% Extract average spectra for each annotation type from each file...
for n = 1:numF
    
    % Load each file
    tmp = open([fold file{n}]);
    
    % Extract all annotated regions
    [mask2,~,~] = desiAnnotationExtract(tmp.dpn);
    mask2(mask2 == 15) = 0;
    
    fx = mask2 > 0;
    
    % Determine average spectrum
    sp = tmp.dpn.d1.sp;
    sz = size(sp);
    sp = reshape(sp,[sz(1)*sz(2) sz(3)]);
    
    mz = tmp.dpn.d1.mz';
    tmp = nanmean(sp(fx,:),1)';
    fx = tmp == 0;
        
    
    avgs{n,1} = [mz(~fx) tmp(~fx)];
    avgs{n,2} = file{n};
    
end

end

