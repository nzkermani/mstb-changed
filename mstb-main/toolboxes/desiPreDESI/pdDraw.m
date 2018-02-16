function [fig] = pdDraw
%pdDraw - draw the window for preDESI

bgCol = 'white';

f0 = findobj('Tag','preDESI');
delete(f0);

fig.fig = figure('Units','normalized',...
    'Position',[0.1 0.25 0.5 0.5],...
    'Toolbar','none',...
    'MenuBar','none',...
    'Tag','preDESI',...
    'Color',bgCol);

% File selection button
fig.fSelect = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.01 0.95 0.08 0.049],...
    'Style','pushbutton',...
    'String','File...',...
    'FontSize',18);

% File display box...
fig.fText = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.1 0.95 0.88 0.045],...
    'Style','text',...
    'String','',...
    'FontSize',18,...
    'HorizontalAlignment','left',...
    'BackgroundColor',bgCol);

% Panel for the ion extraction parts...
fig.pan1 = uipanel('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.01 0.01 0.18 0.88],...
    'BackgroundColor','white');

% Panel for the processing parts...
fig.pan2 = uipanel('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.71 0.01 0.28 0.88],...
    'BackgroundColor','white');

% Axes
fig.ax(1) = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.22 0.05 0.46 0.8]);

fig.ax(2) = imagesc(rand(100,100));
set(fig.ax(1),'XTick',[],...
    'YTick',[],...
    'TickLength',[0 0]);
colormap(gray)


% Draw bits in each panel
[fig] = panel1(fig,bgCol);
[fig] = panel2(fig,bgCol);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = panel1(fig,bgCol)

% So now we add the parts for the ion extraction panel
uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0 0.925 1 0.075],...
    'Style','text',...
    'String','I o n s',...
    'FontSize',30,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'BackgroundColor',bgCol);

% Button for pos / neg
fig.polarity = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0 0.8 0.25 0.1],...
    'Style','pushbutton',...
    'String','NEG',...
    'FontSize',14);

% Box for the ions
fig.ionList = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.5 0.6 0.49 0.3],...
    'Style','edit',...
    'String','Ions',...
    'FontSize',14,...
    'Min',1,...
    'Max',10);

% Box for ppm tolerance
uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.01 0.6 0.23 0.045],...
    'Style','text',...
    'String','� ppm',...
    'FontSize',14,...
    'HorizontalAlignment','center',...
    'BackgroundColor',bgCol);

fig.ionPPM = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.25 0.6 0.24 0.05],...
    'Style','edit',...
    'String','5',...
    'FontSize',14);

% Button to do the extraction
fig.extract = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.01 0.5 0.98 0.08],...
    'Style','pushbutton',...
    'String','Extract',...
    'FontSize',20,...
    'FontWeight','bold');

% Now how about visualising the image once it has been extracted?
uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0 0.40 1 0.075],...
    'Style','text',...
    'String','V i s u a l i s e ',...
    'FontSize',30,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'BackgroundColor',bgCol);

% Another box for the ions
fig.ionDisplay = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.5 0.075 0.49 0.3],...
    'Style','listbox',...
    'String','Ions',...
    'FontSize',14,...
    'Min',1,...
    'Max',10);


% Boxes for min/max intensities
uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.01 0.125 0.23 0.045],...
    'Style','text',...
    'String','Min',...
    'FontSize',14,...
    'HorizontalAlignment','center',...
    'BackgroundColor',bgCol);

fig.minInt = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.25 0.125 0.24 0.05],...
    'Style','edit',...
    'String','0',...
    'FontSize',14);

uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.01 0.075 0.23 0.045],...
    'Style','text',...
    'String','Max',...
    'FontSize',14,...
    'HorizontalAlignment','center',...
    'BackgroundColor',bgCol);

fig.maxInt = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.25 0.075 0.24 0.05],...
    'Style','edit',...
    'String','1',...
    'FontSize',14);

fig.doLog = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.01 0.2 0.48 0.05],...
    'Style','checkbox',...
    'String','Log?',...
    'Value',1,...
    'FontSize',14,...
    'BackgroundColor',bgCol);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = panel2(fig,bgCol)

% Segmentation
uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0 0.925 1 0.075],...
    'Style','text',...
    'String','S e g m e n t s',...
    'FontSize',30,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'BackgroundColor',bgCol);

% Number of sections
uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.01 0.85 0.31 0.045],...
    'Style','text',...
    'String','Sections = ',...
    'FontSize',14,...
    'HorizontalAlignment','center',...
    'BackgroundColor',bgCol);

fig.numSect = uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.34 0.85 0.31 0.05],...
    'Style','edit',...
    'String','4',...
    'FontSize',14);

fig.doSegment = uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.67 0.85 0.31 0.05],...
    'Style','pushbutton',...
    'String','Detect',...
    'FontSize',14);

% % Heading for sections...
% uicontrol('Parent',fig.pan2,...
%     'Units','normalized',...
%     'Position',[0 0.725 1 0.075],...
%     'Style','text',...
%     'String','S e g m e n t s',...
%     'FontSize',30,...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'BackgroundColor',bgCol);

% Empty table...
fig.segTable = uitable('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.01 0.525 0.98 0.3],...
    'ColumnName',{'?','Colour','ID'},...
    'ColumnEditable',[true false true],...
    'FontSize',16,...
    'ColumnWidth',{30 75 200});

% Delete segment
fig.segDelete = uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.01 0.475 0.48 0.05],...
    'Style','pushbutton',...
    'String','Delete Segment',...
    'FontSize',14);

fig.segAdd = uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.51 0.475 0.48 0.05],...
    'Style','pushbutton',...
    'String','Add Segment',...
    'FontSize',14);


% Processing
uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0 0.15 1 0.075],...
    'Style','text',...
    'String','P r o c e s s',...
    'FontSize',30,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'BackgroundColor',bgCol);

uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.01 0.06 0.48 0.045],...
    'Style','text',...
    'String','Processing method: ',...
    'FontSize',14,...
    'HorizontalAlignment','center',...
    'BackgroundColor',bgCol);

fig.procOpt = uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.51 0.06 0.48 0.05],...
    'Style','popupmenu',...
    'String',{'Separate';'Individual'},...
    'Value',1,...
    'FontSize',14);

fig.process = uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.01 0.01 0.98 0.05],...
    'Style','pushbutton',...
    'String','Begin Processing...',...
    'FontSize',14);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%