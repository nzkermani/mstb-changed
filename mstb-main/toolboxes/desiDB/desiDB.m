function desiDB
%desiDB - this is the new and rejuvenated functino for launching a DESI
%toolbox database (for the new toolbox!)

% Let's define the default path here
if ismac
    %defP = '/Volumes/Data/Data/';
    %defP = '/Users/jmckenzi/Dropbox/Imperial/Projects/Multiple Samples/Colo/';%DB/MSA/';
    %defP = '/Users/jmckenzi/Dropbox/Imperial/Projects/Metaspace/Engine Dump/MTSP-Merge/';
    defP = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/STATS/';
else
    defP = '\\zoltan1.med.ic.ac.uk\Data\Data\Colorectal\ColorectalDummy\';
end

% Check that it exists. If not, use current folder
if ~exist(defP,'dir')
    defP = [pwd filesep];
end

% Now run the function to determine the files that will go in the database
[allF] = dbFileFind(defP);

% Get the file information for the table
[tabDat] = dbFileInfo(allF);

% Here we can draw the database
[fig] = dbDraw;
set(fig.fig,'Name',defP);

% Draw the toolbar with the buttons
[fig] = dbToolbar(fig);

% Set the table data
dbTableUpdate(fig,tabDat);

% Now update the callbacks of the buttons and stuff
dbToolbarCallback([],[],fig,defP);

% Here we turn on the options button, and then show the uiTree. This uiTree
% will need to be updated quite significantly
%set(fig.proc,'State','on'); 
%dbUItree(fig.proc,[],fig,[]);

return

% This is the button for the selection of all or none of the files
set(fig.selAll,'Callback', {@dbUIselAll,fig,tabData2});

% Here we turn on the options button, and then show the uiTree. This uiTree
% will need to be updated quite significantly
set(fig.nb12,'State','on'); 
dbUItree(fig.nb12,[],fig,files);

% This is the button that allows the tree to be saved
set(fig.tr.def,'Callback',{@uiTreeSaveDef,fig},'Enable','on');
set(fig.tr.reset,'Callback',{@uiTreeResetDef,fig},'Enable','on');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
