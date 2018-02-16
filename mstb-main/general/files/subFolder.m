function [sf] = subFolder(path)
%previousFolder - return the path of the parent folder

if ~strcmp(path(end),filesep)
    path = [path filesep];
end

% Find slashes
sl = strfind(path,filesep);

% New path
sf = path(sl(end-1)+1:sl(end)-1);

end

