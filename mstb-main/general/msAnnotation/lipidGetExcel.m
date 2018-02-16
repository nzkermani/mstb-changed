function [new] = lipidGetExcel
%lipidGetExcel - get the various lipid files from the Z drive.

% Determine the files in the folder
locn = '/Volumes/Data/Group/Projects/data bases lipid maps/lipidmaps/';
[hits] = fileFinderAll(locn,'xlsx');
hits(1,:) = [];
numF = size(hits,1);

% Structure in which to save everything
lm = struct('LMID',[],'Name',[],...
    'IUPAC',[],'Formula',[],...
    'Mass',[],'Class1',[],...
    'Class2',[],'Class3',[],...
    'Source',[]);

% Loop through each one...
for n = 1:numF
    
    file = [hits{n,1} filesep hits{n,2}];
    if ~exist(file,'file')
        continue;
    end
    
    % Read from Excel
    [~,~,c] = xlsread([hits{n,1} filesep hits{n,2}]);
    
    % Determine numeric indices
    idx = cellfun(@isnumeric,c(:,5));
    if sum(idx) == 0
        disp(hits{n,2});
    end
    
    % Splat the things into the structure
    lm(n).LMID      = c(idx,1);
    lm(n).Name      = c(idx,2);
    lm(n).IUPAC     = c(idx,3);
    lm(n).Formula   = c(idx,4);
    lm(n).Mass      = cell2mat(c(idx,5));
    lm(n).Class1    = c(idx,6);
    lm(n).Class2    = c(idx,7);
    lm(n).Class3    = c(idx,8);
    lm(n).Source    = repmat(hits(n,2),[sum(idx) 1]);
    
end

% Concatenate all of them...
fn = fieldnames(lm);
for n = 1:numel(fn)    
    new.(fn{n}) = vertcat(lm.(fn{n}));   
end

end

