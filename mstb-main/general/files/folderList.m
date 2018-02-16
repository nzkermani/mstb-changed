function [ op ] = folderList(path)
%folderList - return a list of folders in the current one - no subfolders
%are selected...

p = dir(path);

% Names
n = {p.name}';

% Which ones are folders?
isDir = vertcat(p.isdir);

% Find any which contain  a dot .
isDot = ~cellfun(@isempty,strfind(n,'.'));

% Which remain?
fx = isDir & ~isDot;

op = n(fx);

end

