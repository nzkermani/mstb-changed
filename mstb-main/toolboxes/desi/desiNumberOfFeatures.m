function [qty,subF] = desiNumberOfFeatures(fold)
% Determine the number of features in a series of files

if isempty(fold)
    fold = '/Volumes/JSM/DB/PBC/PBC Liver/';
end

% Get all mat files
allF = fileFinderAll(fold,'mat',true);
numF = size(allF,1);

% Storage...
qty = zeros(numF,1);

% Extra info?
subF = cell(numF,1);

for n = 1:numF
    
    % Loop through
    tmp = open([allF{n,1} filesep allF{n,2}]);
    
    if ~isfield(tmp,'dpn')
        continue;
    end
    
    qty(n,1) = numel(tmp.dpn.d1.mz);
    
    % What about the subfolder for classification purposes?
    [~,subF{n,1}] = previousFolder(allF{n,1});
    
end



end