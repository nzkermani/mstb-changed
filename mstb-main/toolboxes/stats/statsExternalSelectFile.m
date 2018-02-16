function statsExternalSelectFile(src,~,~,window)
% statsExternalSelectFile - launch dialog to select a file containing the
% external test set samples

% Define the default location
defP = '/Users/jmckenzi/Dropbox/Imperial/Projects/REIMS Laser/Test Set/';
if ~exist(defP,'dir')
    defP = pwd;
end

% Ask user for file...
[fNam,pNam] = uigetfile({'*.mat'},'Select',defP);
if isnumeric(fNam)
    set(window.file,'String','');
    set(src,'UserData',[]);
    set(window.classCompare,'Enable','off','String','N/A','Value',1);
    return
end

% Open the file
tmp = open([pNam fNam]);

% Simple check to ensure that this is a correct looking file...
if isfield(tmp,'res')
    set(window.file,'String',fNam);
else
    set(window.file,'String','');
    set(src,'UserData',[]);
    set(window.classCompare,'Enable','off','String','N/A','Value',1);
    return
end

% Need to determine the meta data fields against which we will compare the
% classification performance
fn = fieldnames(tmp.res.meta);
fi = find(strcmpi(fn,'histid'));
if ~isempty(fi)
    val = fi;
else
    val = 1;
end
set(window.classCompare,'String',fn,...
    'Enable','on',...
    'Value',val);

% Now we can proceed...
test.path = pNam;
test.file = fNam;
set(src,'UserData',test);
clear tmp

end