function [ op ] = prgImport(fn)
%prgImport - read in CSV file from Progenesis / CPC

if nargin == 0
    fn = ['/Users/jmckenzi/Dropbox/Imperial/Projects/JLA/'...
        'JamesAlexander AT HPOS_combinedData.csv'];
end

% Open the file
fid = fopen(fn,'r');

% Read each line
tl = fgetl(fid);
i = 1;
isHeader = true;

% How big should we make the table? Expected number of observations?
numO = 1000;

% Other counters
numH = 0;

while ischar(tl)
    
    % Strip according to commas
    [heads,~] = regexp(tl,',','split');
    
    % Examine first entry of heads, and if 0 is means that this is the
    % first observation record
    if isHeader && strcmp(heads{1},'0')
        isHeader = false;
    end
    
    % If this is the first row, then we need to determine the headings and
    % separate them from the variables which come after, ideally the first
    % being marked as '0'
    if i == 1
        
        % Find the column of the first variable...
        firstV = find(strcmp(heads,'0'));
        numV = numel(heads) - firstV + 1;
        
        % Create a cell for storing all of the information...
        meta = cell(numO,firstV-1);
        
        % Create somewhere to store variable information
        varInfo = cell(numV,10);
        
    end
    
    % Store the information...
    meta(i,:) = heads(1,1:firstV-1);
    
    % Store variable information
    if isHeader && i > 1
        numH = numH + 1;
        varInfo(:,numH) = heads(1,firstV:end)';        
    end
    
        
    % Finally increase the counter
    i = i + 1;
    disp(int2str(i));
    
    % Get the next line ready to go again
    tl = fgetl(fid);
    
end

fclose(fid);

% Trim the various parts
meta = meta(1:i-1,:);
varInfo = varInfo(:,1:numH);

% Determine the first sample...
fs = find(strcmp(meta(:,1),'0'));
varHeads = meta(2:fs-1,1);
meta(2:fs-1,:) = [];

% Can we csvread the file using the appropriate coordinates?
dataCoord = [numH+1 firstV-1];
op.data = csvread(fn,dataCoord(1),dataCoord(2));

op.varH = varHeads;
op.var = varInfo;
op.metaH = meta(1,:)';
op.meta = meta(2:end,:);

% Make it a little tidier
[op] = prgTidy(op);

end

