function [ fl ] = fileFinder(path,type)
%fileFinder - file all files in a folder

% If no path ask for path
if isempty(path)
    path = uigetdir;
end

% Check it exists...
if ~exist(path,'dir')
    error('Folder not found');
end

% Directory contents
allF = dir(path);

numF = size(allF,1);

fl = cell(numF,1);

for n = 1:numF
    
    % Early quit
    if length(allF(n).name) < length(type)
        continue;
    end
    
    if strcmp(allF(n).name(end-length(type)+1:end),type)
        fl(n,1) = {allF(n).name};
    end
    
end

fl = fl(~cellfun(@isempty,fl));
        


end

