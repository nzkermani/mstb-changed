function [ timeStamp ] = watersTimeStamp(raw)
%watersTimeStamp - read in the header file and determine the  acquisition
%date / time from it

% From the 'raw' file path, determine the name of the header file
head = [raw filesep '_HEADER.TXT'];

% What are we looking for?
td = '$$ Acquired Date: ';
tt = '$$ Acquired Time: ';


% Cell of arbitrary length
txt = cell(10,1);
n = 0;

% Create empty lineID variable
lineIDDate = 0;
lineIDTime = 0;

% Read in the file
fid = fopen(head,'r');
while ~feof(fid) && n < size(txt,1)
    n = n + 1;
    txt{n,1} = fgetl(fid);
    
    % Identify the correct line...
    if length(txt{n,1}) >= length(td)
        
        if strcmp(txt{n,1}(1:length(td)),td)
            lineIDDate = n;
        elseif strcmp(txt{n,1}(1:length(tt)),tt)
            lineIDTime = n;
        end
    end
       
end
fclose(fid);

% Format time...
if lineIDDate > 0    
    i1 = txt{lineIDDate,1}(length(td)+1:end);
    i2 = txt{lineIDTime,1}(length(tt)+1:end);
    timeStamp = datenum([i1 '-' i2],'dd-mmm-yyyy-HH:MM:SS');
else
    timeStamp = [];
end

end

