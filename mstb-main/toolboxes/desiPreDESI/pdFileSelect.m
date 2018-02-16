function pdFileSelect(~,~,fig,defP)
%pdFileSelect...

% Ask the user for a file
[file.nam,file.dir,flag] = uigetfile({'*.imzML; *.dat; *.DAT'},...
    'File Select',defP);

% If no file selected, then quit/do nothing
if flag == 0
    fig.fText.String = 'No file selected';
    return
end

% Determine the file extension
file.ext = fileExtension(file.nam);

% Reformat the .dat files to the .raw file name instead
switch lower(file.ext)
    case 'dat'
        sl = strfind(file.dir,filesep);
        sl = sl(end-1);
        file.nam = file.dir(sl+1:end-1);
        file.dir = file.dir(1:sl);
        file.ext = 'raw';
end

% Update the file in the text bar...
fig.fText.String = file.nam;
fig.fSelect.UserData = file;


% Also need to clear the axis of images, patches, boxes etc...
data = guidata(fig.fig);
set(fig.ax(2),'CData',rand(size(get(fig.ax(2),'CData'))));
if isempty(data)
    return
end
if ~isfield(data,'table')
    return
end
if isempty(data.table)
    return
end

delete([data.table{:,1}]);
delete([data.table{:,2}]);

data.table = [];
fig.segTable.Data = [];

guidata(fig.fig,data);

end

