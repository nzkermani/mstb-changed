function [ all,trim ] = fileFoldNames(path)
%fileFoldNames - export the file/folder names from a directory into a
%simple cell. Also need one without interstitial whitespace

ls = dir(path)

% List all files/folders
all = {ls.name}';
trim = cell(size(all));

% Loop through each one
for n = 1:size(all,1)
    
    tmp = all{n};
    ws = isspace(tmp);
    
    trim{n} = tmp(~ws);
    
end

end

