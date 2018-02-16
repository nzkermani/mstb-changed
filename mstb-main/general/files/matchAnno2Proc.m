function [ match ] = matchAnno2Proc
%matchAnno2Proc - for the previously annotated mat files, find a match to
%newly processed files

annoLoc = 'E:\Data\Olivia\Anno\';
procLoc = 'E:\Data\Olivia\Proc\';

% Get list of anno files
annoF = fileFinderAll(annoLoc,'mat',true);
numF = size(annoF,1);

% What about the processed files?
procF = fileFinderAll(procLoc,'mat',true);

% Store match information
match = cell(numF,4);

% Loop through the files
for n = 1:numF
    
    % Get the folder string
    sl = strfind(annoF{n,1},filesep);
    
    fld = annoF{n,1}(sl(end)+1:end);
    
    % Break up according to the _ characters
    usc = strfind(fld,'_');
    numU = numel(usc);
    if numU == 2        
        p1 = fld(1:usc(1)-1);
        p2 = fld(usc(1)+1:usc(2)-1);
        p3 = fld(usc(2)+1:end);
    elseif numU == 1        
        p1 = fld(1:usc(1)-1);
        p2 = '_'; % guaranteed in all files
        p3 = fld(usc(1)+1:end);        
    end
    
    % Do the matching here
    m1 = ~cellfun(@isempty,strfind(procF(:,2),p1));
    m2 = ~cellfun(@isempty,strfind(procF(:,2),p2));
    m3 = ~cellfun(@isempty,strfind(procF(:,2),p3));
    mm = m1 & m2 & m3;
    mm = find(mm);
    
    % If there is a single match, then take it!
    if numel(mm) == 1
        
        match{n,1} = annoF{n,1};
        match{n,2} = annoF{n,2};
        match{n,3} = procF{mm,1};
        match{n,4} = procF{mm,2};
        
    end
    
    %match(1:n,[2 4])
    
end

% Strip out the empty ones
fx = cellfun(@isempty,match(:,1));

match = match(~fx,:);

end

