function changeClasses(origClass,newClass)
%changeClasses - if given a histID vector of group membership, and some new
%class labels, allow the user to define which classes could be merged and
%write out a config file, i.e. new classes

% Need to decide if this is a structure. If so, then a slighthly different
% operation as need to format the whole thing, data as well...
if isstruct(origClass)
    [unqC,~,~] = unique(origClass.histID);
else    
    [unqC,~,~] = unique(origClass);
end

% Draw a window...
[fig] = drawWindow2(newClass,unqC);

% Set the button's callback
set(fig.go,'Callback',{@doMerge,fig,newClass,unqC,origClass});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = drawWindow2(new,orig)
% The new updated version that uses a table and checkboxes rather than
% radiobuttons for the selection of the groups...

bgCol = 'white';

% How many new classes?
numN = size(new,1) + 2;

% And orig ones?
numO = size(orig,1);

f0 = findobj('Tag','metaBuilder');
close(f0);

% First draw the figure
fig.fig = figure('Name','Meta Builder',...
    'Units','normalized',...
    'Position',[0.125 0.01 0.3 0.9],...
    'Menubar','none',...
    'Toolbar','none',...
    'Number','off',...
    'Tag','metaBuilder',...
    'Color',bgCol);

% Add a box to save the analyses results according to the file name
% specified. Provide a default as per the time of the analysis...
%defName = [dbP 'MetaMerge-' datestr(now,'yymmdd-HHMM')];
% fig.savP = uicontrol('Parent',fig.fig,...
%     'Style','pushbutton',...
%     'Units','normalized',...
%     'Position',[0.01 0.81 0.5-0.01 0.03],...
%     'String',dbP.local,...
%     'FontSize',14,...
%     'Callback',{@changePath});
% fig.savN = uicontrol('Parent',fig.fig,...
%     'Style','edit',...
%     'Units','normalized',...
%     'Position',[0.5+0.01 0.81 0.5-0.02 0.03],...
%     'String',[datestr(now,'yymmdd-HHMM')],...
%     'FontSize',14);
% uicontrol('Parent',fig.fig,...
%     'Style','text',...
%     'Units','normalized',...
%     'Position',[0.01 0.84 0.5-0.01 0.03],...
%     'String','Path Name',...
%     'FontSize',16,...
%     'BackgroundColor',bgCol);
% uicontrol('Parent',fig.fig,...
%     'Style','text',...
%     'Units','normalized',...
%     'Position',[0.5+0.01 0.84 0.5-0.02 0.03],...
%     'String','File Name',...
%     'FontSize',16,...
%     'BackgroundColor',bgCol);

% Big button...
fig.go = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.25 0.925 0.5 0.05],...
    'Style','pushbutton',...
    'String','Merge',...
    'FontSize',14);

% Add the table
fig.tab = uitable('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.01 0.01 0.98 0.79]);

% Need to prepare the table data, number of rows equal to the number of
% orig classes, number of columns from new with the orig/delete options
numO = size(orig,1);
numN = size(new,1);

% The initial false boxes...
fls  = cell(1,numN);
edit = zeros(1,numN);
form = cell(1,numN);
cols = cell(1,numN);
for n = 1:numN    
    fls{1,n} = false;
    edit(1,n) = 1;
    form{1,n} = 'logical';
    cols{1,n} = new{n,1};  
end

% Now make the full size rows
for n = 1:numO
    
    % This is the name / n*false, then true for orig and false for delete
    tmp = {orig{n,1} fls{:} false false};
    
    
    if n == 1
        td = tmp;        
        edit = [0 edit 1 1];
        form = {'char' form{:} 'logical' 'logical'};
        cols = {'' cols{:} 'Original' 'Delete'};
        
    else
        td = cat(1,td,tmp);
    end
end
edit = edit == 1;

% Column widths...
colW = cell(1,numel(cols));
for n = 1:numel(cols)
    if n == 1
        colW{1,1} = 300;
    else
        colW{1,n} = 80;
    end
end

% Update the table with the data, editableness, format and column names
set(fig.tab,'Data',td,...
    'ColumnEditable',edit,...
    'ColumnFormat',form,...
    'ColumnName',cols,...
    'FontSize',16,...
    'ColumnWidth',colW);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doMerge(~,~,fig,new,orig,full)

% Manage for struct/cell input of orig
if isstruct(full)
    all = full;
    full = all.histID;
else
    all = [];
end

% Get the composition of each of the new groups...    
numN = size(new,1);
numO = size(orig,1);
chk  = zeros(numO,numN);
td   = get(fig.tab,'Data');
for n = 2:numN+3
    for r = 1:numO        
        % This is for the radio button approach should you wish to
        % reinstate it later on...
        %chk(r,n) = get(fig.rb(r,n),'Value'); 
        
        % This is for the table form of the window
        chk(r,n-1) = td{r,n};
    end
end

% If no class selected, then delete it by default
nc = sum(chk,2) == 0;
chk(nc,end) = 1;


% Create a new empty cell for the new classifications to go in
newFull = cell(size(full));

% Loop through each of the unique annotations
for n = 1:numO
    
    % Indices of original annotations
    fx = strcmp(full,orig{n,1});
    
    % What is it going to be?
    if sum(chk(n,:)) > 1
        disp(['Double selection: ' orig{n,1} ' set to DELETE']);
        class = {'DELETE'};
    
    elseif chk(n,end) == 1
        class = {'DELETE'};
        
    elseif chk(n,end-1) == 1
        class = orig(n,1);
        
    else
        class = new(chk(n,1:end-2)' == 1);
    end
        
    % Update the class
    if ~isempty(class)
        newFull(fx,1) = class;
    end
    
end

% Now we need to format the full structure output if it was provided
if ~isempty(all)
    
    % These are the observations to keep
    fx = ~strcmp(newFull,'DELETE');
    
    % What is the original size?    
    szO = size(all.X)
    
    % What about the other fields in the structure?
    fn = fieldnames(all)
    numF = numel(fn);
    for n = 1:numF
        
        sz = size(all.(fn{n}));
        
        if sz(1) == szO(1) && sz(2) == 1
            
            tmp = all.(fn{n});            
            all.(fn{n}) = tmp(fx,:);
            
        end
    end
    
    % FInally trim the data matrix
    all.X = all.X(fx,:);
    
    % And all the new class
    all.new = newFull(fx,:);
    
    % Return this to the workspace
    assignin('base','newX',all);
    
else

    % Now that we have everything, then just assign the cell back to the
    % workspace
    assignin('base','reclass',newFull);

end

% Finally close this to make it go away
close(fig.fig);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
