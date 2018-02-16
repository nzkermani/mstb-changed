function [ file ] = dat2raw(file)
%dat2raw - convert a dat file type to the raw filetype instead

% Is it a structure?
if ~isstruct(file)
    return
end

% Continue only if this is correct
if ~strcmpi(file.ext,'dat')
    return
end

% Find all file separators
sl = strfind(file.dir,filesep);

% We only care about the last one
sl = sl(end-1);

% Change to incorporate the root and .raw as the file name
file.nam = file.dir(sl+1:end-1);
file.dir = file.dir(1:sl);
file.ext = 'raw';

end

