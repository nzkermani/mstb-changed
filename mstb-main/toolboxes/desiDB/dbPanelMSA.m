function dbPanelMSA(~,~,fig,defP)
%dbPanelMSA - for multiple sample analysis

% Change the default path
defP = get(fig.fig,'Name');

% Draw the panel's bits and bobs
[pan] = panelBits(fig,defP);

% Add callback to change path
set(pan.path,'Callback',{@changePath});

% Add the callback to determine the options from the panel - the function
% should get the options, and then launch into another function that
% actually does the MSA
set(pan.proc,'Callback',{@dbOptionsMSA,fig,pan});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pan] = panelBits(fig,defP)
% Window with the little options for manipulating the figure...

% This is where we draw everything
parent = fig.panel;

fS = 16;

% Heading
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.95 1 0.05],...
    'Style','text',...
    'String','Multiple sample analysis',...
    'FontSize',24,...
    'BackgroundColor',[1 1 1]);

% m/z tolerance (for histogram)
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.8 0.25 0.1],...
    'Style','text',...
    'String','Resolution / Da',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
pan.mzTol = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.275 0.8 + 0.058 0.15 0.04],...
    'Style','edit',...
    'String','0.01',...
    'FontSize',fS);

% ppm tolerance (for peak matching)
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.75 0.25 0.1],...
    'Style','text',...
    'String','Tolerance / ppm',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
pan.ppmTol = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.275 0.75 + 0.058 0.15 0.04],...
    'Style','edit',...
    'String','20',...
    'FontSize',fS);

% mz Low
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.7 0.25 0.1],...
    'Style','text',...
    'String','mz Low',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
pan.mzLow = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.275 0.7 + 0.058 0.15 0.04],...
    'Style','edit',...
    'String','600',...
    'FontSize',fS);

% mz High
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.65 0.25 0.1],...
    'Style','text',...
    'String','mz High',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
pan.mzHigh = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.275 0.65 + 0.058 0.15 0.04],...
    'Style','edit',...
    'String','1000',...
    'FontSize',fS);

% Ion mode - important for later on (and annotations too)
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.5 0.8 0.25 0.1],...
    'Style','text',...
    'String','Polarity',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
pan.polarity = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.775 0.8 0.2 0.1],...
    'Style','listbox',...
    'String',{'Negative','Positive','Switched'},...
    'FontSize',fS-2,...
    'Value',1,...
    'Enable','off');

% Pixel selection
if strcmp(getUser,'jmckenzi')
    pixSel = 'on';
else
    pixSel = 'off';
end
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.5 0.65 0.25 0.1],...
    'Style','text',...
    'String','Pixel Selection',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
pan.pixel = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.775 0.65 0.2 0.1],...
    'Style','listbox',...
    'String',{'Annotated','Metaspace'},...
    'FontSize',fS-2,...
    'Value',1,...
    'Enable',pixSel);

% Pixel combo
if strcmp(getUser,'jmckenzi')
    poss = {'All-QC','All','Random 5','TT Mean','TT Median'};
else
    poss = {'All','Random 5','TT Mean','TT Median'};
end
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.5 0.5 0.25 0.1],...
    'Style','text',...
    'String','Pixel Combo',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
pan.pixelCombo = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.775 0.5 0.2 0.1],...
    'Style','listbox',...
    'String',poss,...
    'FontSize',fS-2,...
    'Value',1,...
    'Enable','on');

% Save location
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.35 0.25 0.1],...
    'Style','text',...
    'String','Path',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
pan.path = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.275 0.35 + 0.058 0.7 0.04],...
    'Style','pushbutton',...
    'String',[defP filesep],...
    'FontSize',fS);

% Save name
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.3 0.25 0.1],...
    'Style','text',...
    'String','Filename',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
pan.name = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.275 0.3 + 0.058 0.7 0.04],...
    'Style','edit',...
    'String',['desiDB-MSA-' datestr(now,'yymmdd-HHMMSS')],...
    'FontSize',fS);

% Minimum number of pixels per annotated tissue
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.105 0.25 0.1],...
    'Style','text',...
    'String','Min Pix.',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
pan.minPix = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.275 0.105 + 0.058 0.15 0.04],...
    'Style','edit',...
    'String','0',...
    'FontSize',fS);

% Date low
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.055 0.25 0.1],...
    'Style','text',...
    'String','Date Start',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
pan.dateLow = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.275 0.055 + 0.058 0.15 0.04],...
    'Style','edit',...
    'String','010100',...
    'FontSize',fS);

% Date high
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.005 0.25 0.1],...
    'Style','text',...
    'String','Date End',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
pan.dateHigh = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.275 0.005 + 0.058 0.15 0.04],...
    'Style','edit',...
    'String',datestr(now,'ddmmyy'),...
    'FontSize',fS);



% Process button
pan.proc = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.81 0.01 0.18 0.1],...
    'Style','pushbutton',...
    'String','Go',...
    'FontSize',fS);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changePath(src,event)
% Change the path of the folder

% Read the current directory
cud = get(src,'String');

% Ask the user for a new one
ned = uigetdir(cud);

% Set, if not empty
if ~isnumeric(ned)
    set(src,'String',[ned filesep]);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
