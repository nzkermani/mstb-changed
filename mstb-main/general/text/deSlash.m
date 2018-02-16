function [ path ] = deSlash( path )
%deSlash - essentially takes a path and makes it compatiable with the
%computer, as i'd thought that matlab was capable of sorting out a
%cacophany of forward and back slashes.
if ismac
    sl = strfind(path,'\');
    path(sl) = '/';    
else
    sl = strfind(path,'/');
    path(sl) = '\';
end

end

