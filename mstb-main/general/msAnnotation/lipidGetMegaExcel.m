function [new] = lipidGetMegaExcel
%lipidGetExcel - get the various lipid files from the Z drive.

% Determine the files in the folder
locn = '/Volumes/Data/Group/Projects/data bases lipid maps/classification_lipids.xlsx';

% Structure in which to save everything
lm = struct('Formula',[],...
    'Mass',[],'Class1',[],...
    'Class2',[],'Class3',[],...
    'Source',[]);

% I want to know the sheets from the second to the last...
[~,sheets] = xlsfinfo(locn);
sheets = sheets(2:end)';
numS = numel(sheets);

% Loop through each sheet
for n = 1:numS
    
    % Read from Excel
    [a,b,c] = xlsread(locn,sheets{n});
    
    % Define indices
    idx = 2:size(b,1);
    
    % Determine numeric indices
%     idx = cellfun(@isnumeric,c(:,2));
%     idx2 = ~isnan(cell2mat(c(idx,2)))
%     idx3 = idx(idx2)
%     if sum(idx) == 0
%         disp(sheets{n});
%     end
    
    % Splat the things into the structure
    %lm(n).LMID      = c(idx,1);
    %lm(n).Name      = c(idx,2);
    %lm(n).IUPAC     = c(idx,3);
    lm(n).Formula   = c(idx,1);
    lm(n).Mass      = cell2mat(c(idx,2));
    lm(n).Class1    = c(idx,5);
    lm(n).Class2    = c(idx,4);
    lm(n).Class3    = c(idx,3);
    lm(n).Source    = repmat({sheets{n}},[numel(idx) 1]);
    
end

% Concatenate all of them...
fn = fieldnames(lm);
for n = 1:numel(fn)    
    new.(fn{n}) = vertcat(lm.(fn{n}));   
end

end

