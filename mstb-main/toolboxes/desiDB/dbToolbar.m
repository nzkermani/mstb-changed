function [fig] = dbToolbar(fig)
% Draw the toolbar for the new DESI DB

% Let's add a new toolbar to it instead of the horrible push buttons
fig.tb.main = uitoolbar('Parent',fig.fig,...
    'Tag','dbtb',...
    'Visible','on');

% Need to get the icons for the buttons...
dirIcon = deSlash([pwd filesep '/desi/icons/']);

% Here a button for the loading a new file...
fig.tb.new = uipushtool('Parent',fig.tb.main,...
    'CData',getImage(dirIcon,'folder-3-48'),...
    'TooltipString','Change the folder');

% Load sample into toolbox
fig.tb.load = uipushtool('Parent',fig.tb.main,...
    'CData',getImage(dirIcon,'tool-box-48'),...
    'TooltipString','Load file into toolbox...');

% Separate!
insertSeparator(fig.tb.main);
insertSeparator(fig.tb.main);
insertSeparator(fig.tb.main);

% A button to toggle the uiTree for processing - assuming that we actually
% have a similar layout...
fig.tb.msaMake = uitoggletool(...
    'CData',getImage(dirIcon,'copy-48'),...
    'TooltipString','Display processing options');

% This button for loading up an already processed MSA model
fig.tb.msaInspect = uitoggletool(...
    'CData',fliplr(getImage(dirIcon,'train-5-48')),...
    'TooltipString','Display processing options');

% Separate
insertSeparator(fig.tb.main);
insertSeparator(fig.tb.main);
insertSeparator(fig.tb.main);

% Close other windows
fig.tb.duck = uipushtool('Parent',fig.tb.main,...
    'CData',getImage(dirIcon,'duck-48'),...
    'TooltipString','Quack - close other windows');

% What about enlarging the toolbar?
dpnEnlarge(fig.tb.main);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ico] = getImage(path,name)
% Get the necessary icon.  Need to beautify it a little bit more...

[ico] = importdata([path name '.png']);

[ico] = iconProcess(ico.alpha,ico.cdata);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function insertSeparator(h)

uipushtool('Parent',h,...
    'CData', NaN([50 50 3]),...
    'ClickedCallback',{},...
    'Separator', 'on',...
    'Tooltip','Is this button missing?');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
