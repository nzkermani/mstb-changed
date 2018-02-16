function [ pass ] = xxxHasAnnotations
%xxxHasAnnotations - check if a folder's worth of files have annotations

dir = '/Volumes/Data/Lab Data/Endometrial/Olivia/Olivia DESI Raw data/1 scan per sec/';

% Get all files
allF = fileFinderAll(dir,'mat',true);
numF = size(allF,1);

pass = cell(numF,3);

% Loop
for n = 1:numF
    
    tmp = open([allF{n,1} filesep allF{n,2}]);
    
    pass{n,1} = allF{n,1};
    pass{n,2} = allF{n,2};
    
    if isfield(tmp,'dpn')
        if isfield(tmp.dpn,'anno')
            pass{n,3} = 'Yes';
        else
            pass{n,3} = 'None';
        end
    end
    
    disp([int2str(n) '/' int2str(numF)]);
end


end

