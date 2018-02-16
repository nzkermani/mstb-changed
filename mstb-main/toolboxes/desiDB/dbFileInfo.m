function [tabData1] = dbFileInfo(files)
% Get the table data based on the unq structure

% How many files?
numF = size(files,1);

% Headings
colName = {'', 'Path', 'Name'};
colForm = {'logical', 'char', 'char'};
colEdit = [true false false];
colWdth = {20 200 200};

% Create the empty table
td = cell(numF,3);
for n = 1:numF    
    td(n,:) = {false files{n,1} files{n,2}};    
end

% Prepare for output as a single structure
tabData1.data = td;
tabData1.name = colName;
tabData1.form = colForm;
tabData1.edit = colEdit;
tabData1.wdth = colWdth;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
