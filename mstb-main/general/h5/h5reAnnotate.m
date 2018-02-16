function [  ] = h5reAnnotate( src,event )
%h5reAnnotate - based on the imzMLbatch function, this modified version
%aims to re-annotate a bunch of h5 files with the new meta-data in an Excel
%spreadsheet.

% Close existing
fO = findobj('Tag', 'h5reanno');
close(fO);

% Need to check that there is a username...
try
    uName = dbCheckName;
catch
    uName = 'anon';
end

% Draw a gui to display the choices for file processing...
[fig] = drawGUI(uName);

% Get the default path...
tmp = open('defaults.mat');
dbPath = tmp.defaults.dbPath;
clear tmp;

% Update the callbacks
set(fig.imzml, 'Callback', {@uiFileSelector,0,fig,dbPath,...
    {'*.imzML'}});

% set(fig.opimg, 'Callback', {@uiFileSelector,0,fig,dbGetDefaults('paths.img'),...
%     {'*.png; *.PNG; *.jpg; *.JPG; *.jpeg; *.JPEG; *.tif; *.TIF','All image files'}});

set(fig.metaP, 'Callback', {@uiFileSelector,1,fig,dbPath,...
    {'*.xlsx'; '*.xls'}});

set(fig.metaS, 'Callback', {@uiMetaSheet,fig});
set(fig.proc, 'Callback', {@doProc,fig});



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [defT] = getTextMessage(opt)
% Get some standard messages for the status window...
switch opt
    case 0
        defT = ['Pick H5 location'];
        
    case 1
        defT = [char(10) char(10)  char(10) ...
            'Now select the '...
            'location of optical images.  Suitable '...
            'formats are png, jpg, jpeg, tif.'];
        
    case 2
        defT = [char(10) char(10) char(10)...            
            'And now the Excel spreadsheet'];
        
    case 3
        defT = [char(10) char(10) char(10) ...
            'Which sheet in the workbook contains the metadata?'];
        
    case 5
        defT = [char(10) char(10) char(10) ...
            'Before pressing the process button, ensure that the ' ...
            'entries match.  Ticked files below will be processed.'];
        
    case 6
        defT = [char(10) char(10) char(10) ...
            'Processing now...' char(10) datestr(now,'HH:MM:SS')];

    case 7
        defT = [char(10) char(10) char(10) ...
            'Finished processing at' char(10) datestr(now,'HH:MM:SS')];
        
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = drawGUI(name)

% Draw the main window
fig.fig = figure('Name', ['<' name '> H5 re-annotation'], 'Tag', 'h5reanno',...
    'Units', 'normalized', ...
    'Position', [0.55 0.15 0.25 0.5],...
    'Color', 'white',...
    'Menubar','none');

% Relative dimensions of the elements
eN = 12;
eH = 1/eN;
%eW = 0.8;

% Add the elements
fig.imzml = uicontrol('Parent', fig.fig, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized',...
    'Position', [0.6 1-(2*eH) 0.4 eH*2], ...
    'String', 'H5 Path...',...
    'Tag','imzML');

fig.metaP = uicontrol('Parent', fig.fig, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized',...
    'Position', [0.6 1-(3*eH) 0.4 eH], ...
    'String', 'Meta Path...',...
    'Enable','off',...
    'Tag','metaP');

fig.metaS = uicontrol('Parent', fig.fig, ...
    'Style', 'listbox', ...
    'Units', 'normalized',...
    'Position', [0.6 1-(5*eH) 0.4 eH*2], ...
    'String', {'...'; 'Pick an excel'; 'worksheet that'; 'contains the metadata'}, ...
    'Enable', 'off',...
    'Tag','metaS');

fig.proc = uicontrol('Parent', fig.fig, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized',...
    'Position', [0 1-(5*eH) 0.6 eH], ...
    'String', 'Process...',...
    'Enable', 'off',...
    'UserData', 0,...
    'Tag','proc');

fig.text = uicontrol('Parent', fig.fig, ...
    'Style', 'text', ...
    'Units', 'normalized',...
    'Position', [0 1-(4*eH) 0.6 eH*4], ...
    'String', getTextMessage(0),...
    'Enable', 'on',...
    'FontSize', 14,...
    'HorizontalAlignment','left');

fig.tab = uitable('Parent', fig.fig, ...
    'Units', 'normalized',...
    'Position', [0 0 1 7*eH]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uiFileSelector(src,event,flag,fig,defP,fType) %#ok<INUSL>
% Generic function to launch a ui get file thing

if flag == 1
    % This is for a single file
    [p.fN,p.fP,~] = uigetfile(fType, 'Pick a file', defP);
else
    % This is for a folder
    [p.fP] = uigetdir(defP, 'Pick a folder');
end
    
% Check that something was pressed...
if length(p.fP) < 2
    disp('No file selected');
    set(src,'String','< >');
    tmp = get(fig.proc, 'UserData');
    set(fig.proc, 'UserData', max([tmp-1 0]));

    return
end

% Update the name of the button to reflect the new file name
set(src, 'String',   p.fP);
set(src, 'UserData', p);

if flag == 1
    set(src,'String',[p.fP p.fN]);
    
    % Update excel sheet information
    [~,sheets] = xlsfinfo([p.fP p.fN]);

    % Enable the choice
    set(fig.metaS, 'String', sheets, 'Value', 1, 'Enable', 'on');    
    
    % Will want to update the table here...    
end

% Enable the next button...
tag = get(src,'Tag');
switch tag
    case 'imzML'
        set(fig.metaP,'Enable','on');
        set(fig.text,'String',getTextMessage(2));
    case 'metaP'
        set(fig.metaS,'Enable','on');
        set(fig.text,'String',getTextMessage(3));
    case 'metaS'
        set(fig.tab,'Enable','on');
        set(fig.text,'String',getTextMessage(4));
    case 'tab'
        set(fig.proc,'Enable','on');
        set(fig.text,'String',getTextMessage(5));
end

% Add one to the userdata for the main go button
% tmp = get(fig.proc, 'UserData');
% set(fig.proc, 'UserData', tmp+1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uiMetaSheet(src,event,fig) %#ok<INUSL>
% A function to get the sheet from the xl spreadsheet provided

% Clear the table first
set(fig.tab, 'Data',[]);

% What has been picked?
tmp1 = get(src, 'Value');
tmp2 = get(src, 'String');

% Need to read in this sheet...
sheet = tmp2{tmp1};
fName = get(fig.metaP, 'String');

% Read it in ...
[d.num,d.txt,d.raw] = xlsread(fName, sheet);
stp = max([size(d.num,1) size(d.txt,1)]);

for n = 2:stp
    try
        metaNames{n-1,1} = char(d.raw(n,1)); %#ok<AGROW>
    catch
        tt = d.raw(n,1);
        tt = cell2mat(tt);
        metaNames{n-1,1} = int2str(tt); %#ok<AGROW>
    end
    tmp = metaNames{n-1,1};
    if strcmp(tmp(end-2:end),'.h5')
        metaNames{n-1,1} = tmp(1:end-3);
    end
end

set(fig.metaS, 'UserData', metaNames, 'Enable', 'on');
set(fig.proc,'Enable','on');
set(fig.text,'String',getTextMessage(5));

uiTablePop(fig);

% Update the text string...
set(fig.text,'String',getTextMessage(5))

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uiTablePop(fig)
% This function runs when all three things have been done

%charCmp = 3;

% Generate an iterative path list
locn = get(fig.imzml, 'String');
if ismac
    fL = ['x:' genpath(locn)];
    seps = strfind(fL,':');
else
    fL = ['x;' genpath(locn)];
    seps = strfind(fL,';');
end

% Somewhere to store the data...
numF = 0;
inf     = struct('h5',[],'meta',[],'good',[],'path',[],'mNum',[]);
inf2    = struct('h5',[],'meta',[],'good',[],'path',[],'mNum',[]);

% Find the separations...

for n = 1:numel(seps)-1
    i1 = seps(n)  +1;
    i2 = seps(n+1)-1;
    path = [fL(i1:i2) filesep];
    
    % Now get the contents of this directory
    dire = dir(path);
    
    % Run through each possible entry...
    for r = 1:size(dire,1)
        
        % Decide if it is an imzML file...
        if ~dire(r).isdir && ~strcmp(dire(r).name(1),'.') ...
                && strcmp(dire(r).name(end-1:end),'h5')
            h5Name = dire(r).name(1:end-3);
                
            % Set the goodness to 0, only gets to be 1 when we find a
            % metadata entry
            numF = numF + 1;
            inf(numF).h5 = h5Name;
            inf(numF).good = 0;
            inf(numF).path = path;

        end
    end    
end

% Get the name of the entries in the metadata spreadsheet...
meta = get(fig.metaS, 'UserData');

% Now can run through the inf structure and match the entries there
% hopefully to one in the meta structure
for n = 1:numF    
    cmp = strcmp(meta, inf(n).h5);
    if sum(cmp) == 1
        fx = find(cmp == 1);        
        inf(n).meta = meta{fx,1};
        inf(n).mNum = fx;
        inf(n).good = 1;        
    end    
end
    
% Remove empty
nn = 0;
for n = 1:size(inf,2)    
    if inf(n).good == 1
        nn = nn + 1;
        inf2(nn).h5 = inf(n).h5;
        inf2(nn).path = inf(n).path;
        inf2(nn).good = inf(n).good;
        inf2(nn).mNum = inf(n).mNum;
        inf2(nn).meta = inf(n).meta;
    end
end

if nn == 0
    msgTxt = ['No files have been matched. Ensure that everything has '...
        'the exact name and that the images are in the same folder as '...
        'the imzML files.  Also check that you''ve picked the right '...
        'worksheet within the spreadsheet.'];
    msgbox(msgTxt,'Woops! This is embarrassing.','help','modal');
    return
else
    inf = inf2;
end

% Now that we have found our matches...put into the table?
numG = 0;
for n = 1:nn    
    if inf(n).good == 1 && inf(n).mNum ~= 0
        numG = numG + 1;
        tmp = {true inf(n).h5 inf(n).meta};
        
        if numG == 1;
            tabData = tmp;
        else
            tabData = cat(1,tabData,tmp);
        end
    end
end

tabDat.data = tabData;
tabDat.colForm = {'logical', 'char', 'char'};
tabDat.colEdit = [true false false];
tabDat.colWdth = {20 195 195};
tabDat.colName = {'x', 'imzML', 'Meta'};

set(fig.tab, 'Data', tabDat.data,...
    'ColumnName', tabDat.colName, ...
    'ColumnFormat', tabDat.colForm, ...
    'ColumnEditable', tabDat.colEdit, ...
    'ColumnWidth', tabDat.colWdth, ...
    'UserData', inf);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doProc(src,event,fig) %#ok<INUSL>
% Call this function to do the processing and then launch the image toolbox
% and then look to encouage the stuff to be saved...

% Get the status of checked boxes
td = get(fig.tab, 'Data')';
td = cell2mat(td(1,:));
sel = find(td == 1);
numF = numel(sel);

% Quit if none selected...
if numF == 0
    msgbox('No files were selected.','Woops!', 'help','modal');
    return
end

% Update the text field
set(fig.text,'String',getTextMessage(6));


% Useful data from the table
inf = get(fig.tab, 'UserData');

% Now for each file that has been selected...
for n = 1:numel(sel)   
    i = sel(n);
    % Now run the processing file...
    if inf(i).good == 1 && inf(i).mNum ~= 0
        doReAnno(fig,inf(i));
    end
end

% Final update
msg = {get(fig.text,'String') getTextMessage(7)};
set(fig.text,'String',char(msg));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doReAnno(fig,inf)
% This is called by another function that provides all the useful info

% Read the specific meta information...
xFile  = get(fig.metaP, 'String');
xSheet = get(fig.metaS, 'String');
tmp3   = get(fig.metaS, 'Value');
xSheet = xSheet{tmp3};

% Let's read in the meta data that we wish to place into the files... First
% the headers
[~,heads,~] = xlsread(xFile, xSheet, '1:1');
numH = size(heads,2);
tmp1 = inf.mNum;

% This is the actual data for this file...
[d.num,d.txt,d.raw] = xlsread(xFile, xSheet, [int2str(tmp1+1) ':' int2str(tmp1+1)]);
h = dbFormatXLmeta(numH,heads,d);

% Delete the existing metadata atrributes... Potentiall troublesome...
h5delMeta([inf.path inf.h5 '.h5']);

% Now that we have the metadata and the file name, we will start to perform
% the overwriting...
numM = size(h,2);
for n = 1:numM    
    h5MetaModify([inf.path inf.h5 '.h5'],{h(n).h},h(n).v);
end    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
