function [lm] = lipidMaps
%lipidMaps - parse the downloaded LipidMaps database
%
% The DB can be exported in text format using the following URL:
% http://www.lipidmaps.org/rest/compound/lm_id/LM/all/download
%
% Further information about using REST is available here:
% http://www.lipidmaps.org/data/lm_rest.php
%
% There may be easier ways to get the data, but it is a good start...

% Which file?
path = '/Users/jmckenzi/Desktop/LMSD.txt';
maxL = 50000;

% Read it in...
fid = fopen(path,'r');

% Get the header names...
tmp = fgetl(fid);
heads = textscan(tmp,'%s\t');
heads = heads{1};
date = heads{1};
heads = heads(2:end);

% Somewhere to store the information
info = cell(maxL,numel(heads));

% Loop
i = 0;
while ~feof(fid)% && i < 10;
    
    i = i + 1;

    % Extract the line and then parse into sections...
    tmp = fgetl(fid);

    % Convert to double to determine the tab spacings
    tab = double(tmp);
    fx = tab == 9;

    % Replace with unlikely symbol
    tab(fx) = double('£');
    
    % Parse
    data = textscan(char(tab),'%s','delimiter','£');
    data = data{1}(2:end)';
    info(i,1:numel(data)) = data;
    
    disp(int2str(i));
    
end

% Close the file
fclose(fid);

% Trim the matrix
info = info(1:i,:);

% Let's put it in a structure instead of a table...
for n = 1:numel(heads)
    lm.(heads{n}) = info(:,n);

    if strcmp(heads{n},'exactmass')
        tmp = str2double(lm.exactmass);
        lm.exactmass = tmp;
    end
end
lm.date = date;

end

