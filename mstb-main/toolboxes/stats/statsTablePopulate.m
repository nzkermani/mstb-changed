function statsTablePopulate(~,~,fig)
%statsTablePopulate - fill the table with information from the meta
%structure...

% Get the guidata...
sts = guidata(fig.fig);
meta = sts.raw.meta;

% What are the field names?
names = fieldnames(meta);
numN = numel(names);

% Create empty table
sz = size(meta.(names{1}),1);
tab = cell(sz,numN);

% Loop through each column
for n = 1:numN
    try
        tab(:,n) = meta.(names{n});
    catch err        
        if strcmp(err.message,'Conversion to cell from double is not possible.');
            for r = 1:sz
                tab{r,n} = meta.(names{n})(r);%repmat({'Numeric'},[sz 1]);
            end
        else
            disp('statsTablePopulate: metadata problems');
            err
        end
    end    
end
    
% Split the table accordingly
rowName = tab(:,1);
colName = ['?'; names(2:end)];
tabDat = tab;
tabDat(:,1) = {true};
cw = repmat({100},[1 numel(colName)]);
cw{1} = 50;
cw{2} = 200;



% Set this to be the actual data...
set(fig.ax.table,'Data',tabDat,...
    'RowName',rowName,...
    'ColumnName',colName,...
    'ColumnEditable',[true false(1,numel(colName))],...
    'ColumnWidth',cw);

end

