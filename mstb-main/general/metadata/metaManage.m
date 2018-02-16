function metaManage(xl)
%metaManage - determine the groups and type of data within an excel
%spreadsheet, and work out a way to present it


% First we need to determine the names of the groups and their unique
% values
[gui.meta] = determineUnique(xl.a,xl.b,xl.c);

% Now create some kind of GUI for displaying, selecting and filtering the
% various fields of metadata
[fig,gui.selected,gui.isFile] = metaDraw(gui.meta);

% Save stuff in guidata to make it dynamic
guidata(fig.fig,gui);

% Set callbacks
set(fig.heads,'Callback',{@boxClick,fig});
set(fig.values,'Callback',{@valueSelect,fig});
set(fig.go,'Callback',{@doExport,fig});
set(fig.search,'Callback',{@changeFolder});
set(fig.savePath,'Callback',{@changeFolder});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [meta] = determineUnique(~,b,c)
% Determine uniqueness

% What are the headings?
heads = b(1,:)';
numH = numel(heads);
unk = 1;

% Maximum number of rows?
numR = size(b,1);

% Loop through...
for n = 1:numH
    
    % Get the values for this header
    vals = c(2:numR,n);
    
    % Try various conversions if it is likely to be numeric
    try
        numb = cell2mat(vals);
        flag = true;
    catch
        numb = [];
        flag = false;
    end
    
    % Check to see if there are any NaNs in the vector and replace these
    % with the work 'Unknown'
    numCheck = cellfun(@isnumeric,vals,'UniformOutput',false);
    tmp = vals;
    tmp(~cell2mat(numCheck)) = {1};
    nanCheck = cellfun(@isnan,tmp,'UniformOutput',false);
    vals(cell2mat(nanCheck)) = {'Unknown'};
    
    % Find any remaining numbers and convert them to text
    if ~flag
        numCheck = cellfun(@isnumeric,vals,'UniformOutput',false);
        idx = find(cell2mat(numCheck));
        for r = 1:numel(idx)
            vals(idx(r),1) = {num2str(vals{idx(r)})};
        end
    end
    
    % Make the header name something useful
    if length(heads{n}) == 0  %#ok<ISMT>
        heads{n} = ['Unknown_' int2str(unk)];
        unk = unk + 1;
    end
    newhead = strSwap(heads{n},' ','_');
    newhead = strDelete(newhead,'.');
    newhead = strSwap(newhead,'%','P');
    newhead = strDelete(newhead,'?');
    newhead = strSwap(newhead,'/','_');
    newhead = strSwap(newhead,'-','_');
    newhead = strSwap(newhead,'&','_');
    
    % Save the results to a structure?
    if flag
        meta.(newhead) = numb;
    else
        meta.(newhead) = vals;
    end
    
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig,selected,isFile] = metaDraw(meta)
% Draw something

fS = 14;

fig.fig = figure('Name','MetaManage',...
    'Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5]);

% List box for the fields
heads = fieldnames(meta);
fig.heads = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'POsition',[0.01 0.2 0.23 0.79],...
    'Style','listbox',...
    'String',heads,...
    'Value',1,...
    'FontSize',fS);

% List box for entries of each meta group
isFile = true(size(meta.(heads{1}),1),numel(heads));
selected = cell(numel(heads),1);
selected{1} = 1:numel(unique(meta.(heads{1})));
fig.values = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'POsition',[0.26 0.2 0.24 0.79],...
    'Style','listbox',...
    'String',unique(meta.(heads{1})),...
    'Value',selected{1},...
    'FontSize',fS,...
    'Min',1,...
    'Max',3,...
    'ListBoxTop',1);

% Another list box to display the actually selected file names...
ff = strcmpi(heads,'file');
fig.files = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'POsition',[0.51 0.2 0.34 0.79],...
    'Style','listbox',...
    'String',meta.(heads{ff}),...
    'FontSize',fS);

fig.quantity = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'POsition',[0.86 0.95 0.14 0.05],...
    'Style','text',...
    'String',int2str(numel(meta.(heads{ff}))),...
    'FontSize',fS);

% What about various lock mass correction options for Ed?
fig.lm(1) = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'POsition',[0.86 0.8 0.14 0.05],...
    'Style','checkbox',...
    'String','LeuEnk',...
    'FontSize',fS);

fig.lm(2) = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'POsition',[0.86 0.75 0.14 0.05],...
    'Style','checkbox',...
    'String','699.5',...
    'FontSize',fS);

% Button to "go"
fig.go = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.9 0.0 0.1 0.2],...
    'Style','pushbutton',...
    'String','Go',...
    'FontSize',fS);

% Buttons for default folder to search for files
path = '/Users/jmckenzi/Dropbox/Imperial/Projects/REIMS Breast iKnife 2/Spectra/';
if ~exist(path,'dir')
    path = pwd;
end
fig.search = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.01 0.12 0.58 0.06],...
    'Style','pushbutton',...
    'String',path,...
    'FontSize',fS-4);

path = '/Users/jmckenzi/Dropbox/Imperial/Projects/REIMS Breast iKnife 2/Models/';
if ~exist(path,'dir')
    path = pwd;
end
fig.savePath = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.01 0.02 0.58 0.06],...
    'Style','pushbutton',...
    'String',path,...
    'FontSize',fS-4);

fig.saveName = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.61 0.02 0.23 0.06],...
    'Style','edit',...
    'String',['Model-' datestr(now,'yymmdd-HHMMSS')],...
    'FontSize',fS-4);



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function boxClick(~,~,fig)
% Update the box when clicked

% Get the guidata
gui = guidata(fig.fig);

% Get the actual clicked value
val = get(fig.heads,'Value');
heads = get(fig.heads,'String');

% Clicked?
click = heads{val};

% Get any previously selected values
prev = gui.selected{val};

% Get the meta values
tmp = gui.meta.(click);
if isempty(prev)
    
    if iscell(tmp)
        prev = 1:size(unique(tmp),1);
    else
        prev = 1:size(unique(tmp,'rows'),1);
    end
end

% Find unique values
if iscell(tmp)
    tmp = unique(tmp);
else
    tmp = unique(gui.meta.(click),'rows');
end

% Update the box
set(fig.values,...
    'String',tmp,...
    'Value',prev);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function valueSelect(src,~,fig)
% What happens when we click

% Get the guidata
gui = guidata(fig.fig);

% Get the current valeus of selection
sels = get(src,'Value');

% Determine which entry this was
val = get(fig.heads,'Value');

% Put into the structure...
gui.selected{val} = sels;

% Save to the guidata
guidata(fig.fig,gui);

% RUn the filtration function
filter2([],[],fig,val);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filter2(src,event,fig,val)
% Hopefully a simpler more dynamic function

% Get the gui data
gui = guidata(fig.fig);

% Get the values from the middle box that are selected
str = get(fig.values,'String');
idx = get(fig.values,'Value');
str = str(idx,:);

% Now we need to compare against the appropriate entries
heads = fieldnames(gui.meta);
comp = gui.meta.(heads{val});
numS = size(str,1);
fx = zeros(size(comp,1),numS);
for n = 1:numS
    
    if iscell(str)
        fx(:,n) = strcmp(comp,str{n});
    elseif ischar(str)
        %fx(:,n) = strcmp(comp,str(n,:));
        ttt = bsxfun(@eq,comp,str(n,:));
        fx(:,n) = sum(ttt,2) == size(ttt,2);
    elseif isnumeric(str)
        fx(:,n) = comp == str(n);
    end
    
end

% Combine fx into single column
fx = sum(fx,2) > 0;

% Add to the matrix
gui.isFile(:,val) = fx;

% Find those which are all true
all = sum(gui.isFile,2) == size(gui.isFile,2);

% We save the selected files to make it easier later on to extract the
% metadata when it comes to exporting the files
gui.allSelected = find(all);

% Update
set(fig.files,'String',gui.meta.File(all),'Value',1);

% How many files?
numF = sum(all);
set(fig.quantity,'String',int2str(numF));

% Update guidata
guidata(fig.fig,gui);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doExport(src,event,fig)
% Find the files in the appropriate folder according to the lock mass
% naming convention, and then prepare the data in format suitable to run
% various models.

% Get the guidata
gui = guidata(fig.fig);

% Get a list of the required files from the box
fList = get(fig.files,'String');

% Determine where we should be looking regarding the file storage - when
% released into the wild, this will need to be dynamic.
path = get(fig.search,'String');
% path = ['/Users/jmckenzi/Dropbox/Imperial/Projects/'...
%     'REIMS Breast iKnife 2/Spectra/'];
pathLEtrue = 'LeuEnk True/';
pathLEfalse = 'LeuEnk False/';

% Now determine which lock mass corrected files we want from the tick boxes
lmLE = get(fig.lm(1),'Value');
lmNA = get(fig.lm(2),'Value');

% Now get all files from the LEUENK LM folder
if lmLE
    [leHits] = fileFinderAll([path pathLEtrue],'txt');
    leHits(1,:) = [];
    leHits(:,3) = reFormatFileName(leHits(:,2));
    
    % Perform file matching
    [matchLE] = fileMatch(fList,leHits);
end

% Now get all files from the nonLEUENK LM folder
if lmNA
    [naHits] = fileFinderAll([path pathLEfalse],'txt');
    naHits(1,:) = [];
    naHits(:,3) = reFormatFileName(naHits(:,2));
    
    % Perform file matching of the non LE LM corrected files
    [matchNA] = fileMatch(fList,naHits);
end

% So now with these files, we would just look to read them in and assemble
% the various bits of metadata as required...
if lmLE && lmNA
    error('Cannot select both sets of data');
elseif lmLE
    allFiles = matchLE;
elseif lmNA
    allFiles = matchNA;
elseif ~lmNA && ~lmLE
    error('Select a lock mass value');
end

% Now we should perhaps metion those files which haven't been found.  So
% show these to the user, and then let them decide if they want to continue
ff = strcmp(allFiles(:,1),'No File Found') | ...
    strcmp(allFiles(:,1),'Multiple Files');

if sum(ff) > 0
    clc;
    disp('**************************************************');
    disp('**************************************************');
    disp('***********SOME FILES WERE NOT FOUND *************');
    disp('**************************************************');
    disp('**************************************************');
    disp(fList(ff));
    disp('**************************************************');
    disp('**************************************************');
    disp('**************************************************');
    disp('**************************************************');
    disp(['Missing files = ' int2str(sum(ff))]);
    choice = questdlg('Files are missing. Continue?',...
        'QUESTION!',...
        'Yes','No','No');
    
    if strcmp(choice,'No')
        return
    else
        allFiles(ff,:) = [];
        selFiles = gui.allSelected;
        selFiles(ff,:) = [];
    end
else
    selFiles = gui.allSelected;
end

if size(allFiles,1) == 0
    disp('No files');
    return
end

% Run a function to import these files and then formulate the metadata as
% appropriate!
extractEverything([],[],fig,allFiles,selFiles);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new] = reFormatFileName(fn)
% Remove some of the extra bits from the file names to make a better match,
% and append with a .raw heading in order to make a perfect match between
% the files

numF = size(fn,1);
new = cell(numF,1);
for n = 1:numF
    tmp = strfind(fn{n},'__');
    new{n,1} = [fn{n}(1:tmp(end)-1) '.raw'];
end    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [info] = fileMatch(files,allFiles)
% Match the files that we have selected to those in the correct folder...

numF = size(files,1);

info = cell(numF,3);

for n = 1:numF
    
    fx = strcmpi(allFiles(:,3),files{n});
    
    if sum(fx) == 0
        info(n,1) = {'No File Found'};
        
    elseif sum(fx) == 1
        info(n,:) = allFiles(fx,:);
    elseif sum(fx) > 1
        info(n,1) = {'Multiple Files'};
    end
    
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function extractEverything(src,event,fig,fileInfo,fileIdx)
% Import the files into memory and extract the appropriate metadata

gui = guidata(fig.fig);

numF = size(fileInfo,1);

% Create a waitbar because why not
wb = waitbar(1/numF,'Preparing');

% Assemble the metadata by trimming out the extraneous stuff
meta = gui.meta;
fld = fieldnames(meta);
numM = numel(fld);
for n = 1:numM
    
    meta.(fld{n}) = gui.meta.(fld{n})(fileIdx,:);
    
end

% Read in the first file...
tmp = dlmread([fileInfo{1,1} filesep fileInfo{1,2}]);

% Extrac the m/z vector and make a matrix of appropriate size
mz = tmp(:,1);

% Matrix
sp = zeros(numF,numel(mz));
sp(1,:) = tmp(:,2)';

% Read in the files...
for n = 2:numF
    
    % Update waitbar
    waitbar(n/numF,wb);
    
    % Read in the files...
    tmp = dlmread([fileInfo{n,1} filesep fileInfo{n,2}]);
    sp(n,:) = tmp(:,2)';
    
end
    

% Simple TIC normalisation for ensuring that the intensities are at least
% comparable
tic = nansum(sp,2);
sp = bsxfun(@rdivide,sp,tic) * 1000;

% Now we need to save everything relevant
savePath = get(fig.savePath,'String');
if ~exist(savePath,'dir')
    mkdir(savePath);
end
saveName = get(fig.saveName,'String');
save([savePath saveName],'mz','sp','meta');

delete(wb);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeFolder(src,~,~)
% Change the folder of a push button

def = get(src,'String');

pt = uigetdir(def,'Select new path');

if isnumeric(pt)
    return
end

% Set the string to be the folder
set(src,'String',[pt filesep]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%