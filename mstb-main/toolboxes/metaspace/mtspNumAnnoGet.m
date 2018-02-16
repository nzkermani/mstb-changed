function [ numAnno,group ] = mtspNumAnnoGet
%mtspNumAnnoGet - determine the number of annotations

% Where?
%locn = '/Volumes/JSM/DB/Metaspace/EngineDumpOLD/MTSP-Merge/Neg/Breast/';
locn = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Merge/Breast/';

% Find all of the mat files in this place
allF = fileFinderAll(locn,'mat',true);
numF = size(allF,1);

% Storage
numAnno = zeros(numF,1);
group = cell(numF,2);

% Loop
for n = 1:numF
    
    % What group is this file?
    [~,b] = previousFolder(allF{n,1});
        
    % Open the file?
    tmp = open([allF{n,1} filesep allF{n,2}]);
    
    if ~isfield(tmp,'dpn')
        continue;
    end
    
    % numAnno?
    numAnno(n,1) = numel(tmp.dpn.mtsp.mz);

    % Save info
    group{n,1} = allF{n,2};
    group{n,2} = b;
    
    disp(int2str(n));
    
end

fx = cellfun(@isempty,group(:,2));

numAnno = numAnno(~fx,:);
group = group(~fx,:);

end

