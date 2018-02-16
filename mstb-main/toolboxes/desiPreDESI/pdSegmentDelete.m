function pdSegmentDelete(~,~,fig)
%pdSegmentDelete - remove the selected segments from the table and image

% Get the guidata
data = guidata(fig.fig);
if isempty(data)
    return
end

% Next get values from the table that are selected
tabDat = fig.segTable.Data;
sel = cell2mat(tabDat(:,1));
if sum(sel) == 0
    return
end

% So how do we go about deleting it?
data.table

% Delete the line and the patch
delete([data.table{sel,1}])
delete([data.table{sel,2}])

% Now delete the entry in the table
data.table = data.table(~sel,:);

% Update the actual table
fig.segTable.Data = data.table(:,3:5);

% Finally update the guidata
guidata(fig.fig,data);

end

