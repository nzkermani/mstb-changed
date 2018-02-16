function dbUpdateTableMSA(~,~,pan,allHist)
%dbUpdateTableMSA - change the values in the table to match those extra
%groups added by the user

% Get the groups in the box
grps = get(pan.groups,'String');

% Get the table data (if it exists!)
td = get(pan.table,'Data');
if isempty(td)
    % Make the table for the first time
    [td,type,colHeads] = generateTableData(allHist,grps);
    set(pan.table,'Data',td,...
        'ColumnEditable',true,...
        'ColumnFormat',type,...
        'ColumnName',colHeads,...
        'RowName',allHist);
    return    
end

% Otherwise we need to just update the number of columns...
if max(size(grps)) == 1
    if strcmp(grps,'')
        set(pan.table,'Data',td(:,end),...
            'ColumnName','Original');
        return
    end
end

% Determine column names
coln = get(pan.table,'ColumnName');
newT = td(:,end);
sz = size(td,1);
for n = numel(grps):-1:1
    
    match = strcmp(coln,grps{n});
    
    if sum(match) == 0
        % Then a new entry
        newT = cat(2,repmat({false},[sz 1]),newT);
    else
        newT = cat(2,td(:,match),newT);
    end
end
set(pan.table,'Data',newT,...
    'ColumnName',cat(1,grps,{'Original'}),...
    'ColumnEditable',true);
    
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [td,type,colHeads] = generateTableData(allHist,grps)
% Generate the table data for the first time

numOrig = numel(allHist);
numNew  = numel(grps);

td = repmat({false},[numOrig numNew+1]);
td(:,end) = {true};
type = repmat({'logical'},[1 numNew+1]);

colHeads = cat(1,grps,{'Original'});

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%