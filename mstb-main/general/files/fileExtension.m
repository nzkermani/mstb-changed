function [ extn ] = fileExtension(fn)
%fileExtension - determine the extension of the filename provided

% Find the (last) full stop
dot = strfind(fn,'.');
if isempty(dot)
    extn = '';
    return
end

dot = dot(end);

extn = fn(dot+1:end);

end

