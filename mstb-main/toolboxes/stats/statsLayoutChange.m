function statsLayoutChange(src,~,fig,datatype)
%statsLayoutChange - switch between the two layout modes

% Get the state of the button
if ~isempty(src)
    value = get(src,'State');

    switch value

        case 'on' % 'detail'
            % ...if on the change to the max detail view
            [locn] = statsLocations('table',datatype);       

        case 'off' % 'simple'
            % ...if off change to the two-grid view
            [locn] = statsLocations('axes',datatype);

    end
else
    % Don't need to do when initialising the thing with a different
    % datatype
    value = 'off';
    [locn] = statsLocations('axes',datatype);
end

% Delete children within the axes
if strcmp(value,'on')
    delete(fig.ax.scatter.Children);
    legend(fig.ax.scatter,'off');
    delete(fig.ax.conf.Children);
    delete(fig.ax.load.Children);
    
    if strcmp(datatype,'lcms')
        delete(fig.ax.spec.Children);
    end
end

% Change the plots!
fig.ax.spec.Visible = locn.spec.vis;
fig.ax.scatter.Visible = locn.scatter.vis;
fig.ax.conf.Visible = locn.conf.vis;
fig.ax.load.Visible = locn.load.vis;

fig.ax.spec.Position = locn.spec.pos;
fig.ax.scatter.Position = locn.scatter.pos;
fig.ax.conf.Position = locn.conf.pos;
fig.ax.load.Position = locn.load.pos;

% Show/hide the metadata table
fig.ax.table.Visible = locn.tab.vis;

% No need to continue if this is the default layout change
if isempty(src)
    return
end

% We may also wish to draw the side menu for this function as well...
if strcmp(value,'on')
    sts = guidata(fig.fig);
    [man] = manipulateWindow(fig,sts);
    
    % Put the callbacks for the various functions
    set(man.newCreate,'Callback',{@createNewMetaGroup,fig,man});
    set(man.groups,'Callback',{@updateUnique,fig,man});
    set(man.unique,'Callback',{@filterUnique,fig,man});
    set(man.setAction,'Callback',{@writeNewInfo,fig,man});
    set(man.miniSearch,'Callback',{@miniSearchCallback,fig,man});
    set(man.importExtra,'Callback',{@statsAddMoreMeta,fig,man});

else
    % Set the values in the table to be true for all observations...
    td = fig.ax.table.Data;
    td(:,1) = {true};
    fig.ax.table.Data = td;
    
end

% Perhaps we need to plot the LC-MS data again...
if strcmp(value,'off') && strcmp(datatype,'lcms')
    statsPlotSpectra([],[],fig,'proc');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [man] = manipulateWindow(fig,sts)
% Window with the little options for manipulating the figure...

% This is where we draw everything
parent = fig.pan2;
f0 = findobj('Parent',parent);
delete(f0);

fS = 14;

% Heading
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.95 1 0.05],...
    'Style','text',...
    'String','Metadata Labelling',...
    'FontSize',24,...
    'BackgroundColor',[1 1 1]);

%%%%


% Subeading
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.90 1 0.05],...
    'Style','text',...
    'String','Select observations',...
    'FontSize',20,...
    'BackgroundColor',[1 1 1]);


% Need to make a window for the meta groups to be listed here...
x = 0.85;
names = fieldnames(sts.raw.meta);
idx = strcmp(names,'hist') | strcmp(names,'histID');
if sum(idx) == 1
    idx = find(idx);
else
    idx = 1;
end
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 x 0.5 0.05],...
    'Style','text',...
    'String','Groups',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
man.groups = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.5 x-0.15 0.475 0.2],...
    'Style','listbox',...
    'String',names,...
    'Value',idx,...
    'FontSize',fS,...
    'Min',1,...
    'Max',1);

% Need to make a window for the meta groups to be listed here...
x = 0.65;
unq = unique(sts.raw.meta.(names{idx}));
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 x 0.5 0.05],...
    'Style','text',...
    'String','Values',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
man.unique = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.5 x-0.15 0.475 0.2],...
    'Style','listbox',...
    'String',unq,...
    'Value',1,...
    'FontSize',fS,...
    'Min',1,...
    'Max',10);

man.miniSearch = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.04 x-0.025 0.4 0.03],...
    'String','',...
    'Style','edit',...
    'FontSize',fS);


% Subeading
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.4 1 0.05],...
    'Style','text',...
    'String','Set new values',...
    'FontSize',20,...
    'BackgroundColor',[1 1 1]);


% What about showing stuff in the other axes?
x = 0.35;
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 x 0.5 0.05],...
    'Style','text',...
    'String','Metadata group',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
man.setGroup = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.5 x-0.15 0.475 0.2],...
    'Style','popupmenu',...
    'String',names,...
    'Value',1,...
    'FontSize',fS);%'Min',1,...'Max',3);

% Set to what value?
x = 0.325;
man.setAction = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 x 0.475 0.03],...
    'Style','pushbutton',...
    'String','Set new value...',...
    'FontSize',fS,...
    'BackgroundColor',[0.9 0.8 0.9]);
man.setText = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.5 x 0.475 0.03],...
    'String','fdsf',...
    'Style','edit',...
    'FontSize',fS);

% Subeading
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.2 1 0.05],...
    'Style','text',...
    'String','Create new group',...
    'FontSize',20,...
    'BackgroundColor',[1 1 1]);

% Need to make a window for the meta groups to be listed here...
x = 0.175;
man.newName = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.5 x 0.475 0.03],...
    'Style','edit',...
    'String','New-Group-Name',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
man.newCreate = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 x 0.475 0.03],...
    'Style','pushbutton',...
    'String','New group...',...
    'FontSize',fS,...
    'BackgroundColor',[0.95 0.95 0.95]);


% Subeading
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.1 1 0.05],...
    'Style','text',...
    'String','Import extra metadata',...
    'FontSize',20,...
    'BackgroundColor',[1 1 1]);

% Need to make a window for the meta groups to be listed here...
x = 0.08;
man.importExtra = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 x 0.475 0.03],...
    'Style','pushbutton',...
    'String','Select csv/xlsx...',...
    'FontSize',fS,...
    'BackgroundColor',[0.95 0.95 0.95]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function createNewMetaGroup(~,~,fig,man)
% Create a new metadata label with the appropriate name. It is a vector of
% blank values to start with, which will be subsequently added to...

% Need to get the name from the box...
name = man.newName.String;

% Format it so that it contains normal characters, no spaces
alp = isstrprop(name,'alpha');
name = name(alp);

% Set the name back in the original box
man.newName.String = name;

% Now we need to actually create it, put it in the structure, and then
% update the table to show it
sts = guidata(fig.fig);
sts.raw.meta.(name) = repmat({'Unknown'},[size(sts.raw.sp,1) 1]);

% Will also copy to the sts.proc.meta.  Need to work out rationale for
% keeping both sets of metadata
warning('Copying metadata to raw & proc');
sts.proc.meta = sts.raw.meta;

% Save to the guidata
guidata(fig.fig,sts);

% Update the table
statsTablePopulate([],[],fig);

% Update the groups box to include the new group
names = fieldnames(sts.raw.meta);
%man.groups.Value = 1;
man.groups.String = names;

man.setGroup.Value = 1;
man.setGroup.String = names;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateUnique(src,~,fig,man)
% Update new entries for each meta group

% Determine the name of the selected group...
val = src.Value;
names = src.String;
name = names{val};

% Determine unique groups
sts = guidata(fig.fig);
unq = unique(sts.raw.meta.(name));

% Update...
man.unique.Value = 1;
man.unique.String = unq;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filterUnique(~,~,fig,man)
% Select the observations in the table that match the selected ones

vals = man.unique.Value;
strs = man.unique.String;

% Determine the group from the other box
grp = man.groups.String(man.groups.Value);

% Get the meta data
sts = guidata(fig.fig);
list = sts.raw.meta.(grp{1});

% Now we need to match values in list to those selected by the strs/vals
% combination
match = repmat({false},size(list));
if iscell(list)
    for n = 1:numel(vals)    
        tmp = strcmp(list,strs{vals(n)});    
        match(tmp) = {true};    
    end
elseif isnumeric(list)
    
    % Convert strs from text to number...then match as usual
    strs = str2num(strs);     %#ok<ST2NM>
    for n = 1:numel(vals)
        tmp = list == strs(vals(n));
        match(tmp) = {true};
    end    
else
    error('Unknown type');
end

% Now set these to be true in the table...
td = fig.ax.table.Data;
td(:,1) = match;
fig.ax.table.Data = td;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeNewInfo(~,~,fig,man)
% Write text from the box into the appropriate part of the table...

% Get the text...
txt = man.setText.String;

% Which group is it going into?
grp = man.setGroup.String{man.setGroup.Value};

% Prevent overwriting histID
if strcmp(grp,'histID') || strcmp(grp,'fileID')
    warning('Cannot overwrite fileID / histID');
    warndlg('Cannot overwrite fileID OR histID - pick/create another group');
    return
end

% Get the current metadata
sts = guidata(fig.fig);

% Into which values is it actually being written (logical vector)
td = fig.ax.table.Data;
fx = cell2mat(td(:,1));

% Extract metadata field
tmp = sts.raw.meta.(grp);

% Write new parts
tmp(fx) = {txt};
sts.raw.meta.(grp) = tmp;

% Set miniSearch to blank
man.miniSearch.String = '';

% Save in guidata
guidata(fig.fig,sts);

% Update the table...
statsTablePopulate([],[],fig);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function miniSearchCallback(~,~,fig,man)
% Perform a miniSearch callback function 

% Get query string
qry = man.miniSearch.String;

% Get the group that is selected...
grp = man.groups.String(man.groups.Value);

% Get the meta data
sts = guidata(fig.fig);
list = lower(sts.raw.meta.(grp{1}));

% Perform a simple exact match
fx = ~cellfun(@isempty,strfind(list,lower(qry)));
fy = cell(size(fx));
for n = 1:size(fy,1)
    if fx(n,1)
        fy(n,1) = {true};
    else
        fy(n,1) = {false};
    end
end

% Now set these to be true in the table...
td = fig.ax.table.Data;
td(:,1) = fy;
fig.ax.table.Data = td;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
