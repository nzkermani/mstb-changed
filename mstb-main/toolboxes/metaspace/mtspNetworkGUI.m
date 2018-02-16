function mtspNetworkGUI(nw)
%mtspNetworkGUI - display a network in a range of ways...

close all;

% What are the group names?
grpNames = fieldnames(nw.gr.Nodes);
fx = strcmp(grpNames,'Properties');
grpNames = grpNames(~fx,:);

% What is the largest value of node connectivity?
maxConn = max(nw.conn(:));

% Draw the gui
[fig] = drawGUI(grpNames,maxConn);

% Save and format guidata
guidata(fig.fig,nw);

% Add callback functions
set(fig.nwDraw,'Callback',{@calculateNetworkXY,fig});

set(fig.marker,'Callback',{@drawNetworkNodes,fig});
set(fig.colour,'Callback',{@drawNetworkNodes,fig});

%set(fig.edgeMin,'Callback',{@drawNetworkEdges,fig});
set(fig.edgeMax,'Callback',{@drawNetworkEdges,fig});

set(fig.transSlider,'Callback',{@changeEdgeTransparency,fig});
set(fig.edgeMinSlider,'Callback',{@changeEdgeMinConnectivity,fig});

% Update the callback for the scroll wheel across the figure
set(fig.activateSlider,'Callback',{@sliderCallback,fig});
set(fig.fig,'WindowScrollWheelFcn',{@sliderCallback,fig});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = drawGUI(grpNames,maxConn)
% Draw the gui

% What about the font size?
fig.fS = 14;

% Figure
fig.fig = figure('Name','NetworkGUI',...
    'Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5]);

% Panel for controls
fig.pan1 = uipanel('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0 0 0.25 1]);


% Panel for axes + axes
fig.pan2 = uipanel('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.25 0 0.75 1],...
    'BackgroundColor','white');
fig.ax = axes('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.1 0.1 0.8 0.8]);
box(fig.ax,'on');
set(fig.ax,'XTick',[],'YTick',[]);

% Colour bar?
fig.cb = colorbar(fig.ax,...
    'Orientation','horizontal',...
    'Location','SouthOutside',...
    'FontSize',fig.fS);

% Controls for the network

% How to scatter the points?
textLabel(fig.pan1,[0.05 0.9 0.4 0.05],'Method',fig.fS);
fig.method = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.55 0.9 0.4 0.05],...
    'Style','popupmenu',...
    'String',{'tSNE';'PCA';'Auto';'Force';'Circle'},...
    'Value',1,...
    'FontSize',fig.fS);

% tSNE options...
textLabel(fig.pan1,[0.05 0.85 0.4 0.05],'tSNE: dims / perp',fig.fS);
fig.tsneDims = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.55 0.85 0.2 0.05],...
    'Style','edit',...
    'String',30,...
    'FontSize',fig.fS);
fig.tsnePerp = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.75 0.85 0.2 0.05],...
    'Style','edit',...
    'String',10,...
    'FontSize',fig.fS);

% Button to draw the network
fig.nwDraw = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.55 0.775 0.4 0.05],...
    'Style','pushbutton',...
    'String','Draw',...
    'FontSize',fig.fS);

% Marker shape group
textLabel(fig.pan1,[0.05 0.65 0.4 0.05],'Group Marker',fig.fS);
fig.marker = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.55 0.65 0.4 0.05],...
    'Style','popupmenu',...
    'String',grpNames,...
    'Value',min([2 numel(grpNames)]),...
    'FontSize',fig.fS);

% Colour group box
textLabel(fig.pan1,[0.05 0.60 0.4 0.05],'Group Colour',fig.fS);
fig.colour = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.55 0.45 0.4 0.2],...
    'Style','listbox',...
    'String',grpNames(3:end,1),...
    'Value',1,...
    'FontSize',fig.fS,...
    'Min',1,...
    'Max',4);

% Things for the line widths and transparency properties...

% Edge width minimum and scaling factor
textLabel(fig.pan1,[0.05 0.35 0.4 0.05],'Max edge width',fig.fS);
fig.edgeMax = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.55 0.35 0.4 0.05],...
    'Style','edit',...
    'String',5,...
    'FontSize',fig.fS);

% Slider for alpha transparency (easy start)
fig.transValue = textLabel(fig.pan1,[0.05 0.25 0.4 0.05],'Trans. = 0.75',fig.fS);
fig.transSlider = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.55 0.25 0.4 0.05],...
    'Style','slider',...
    'Min',0,...
    'Max',1,...
    'Value',0.75,...
    'SliderStep',[0.05 0.1]);

% Slider for minimum number of connections to be displayed... This is
% harder
md = mod(maxConn,5);
if md ~= 0
    maxConn = maxConn + (5-md);
end

% Determine the default value as being the 75% of maxConn to make the
% initial drawing quicker
initConn = maxConn * 0.75;
mdI = mod(initConn,5);
if mdI ~= 0
    initConn = initConn + (5-mdI);
end    
fig.edgeMinValue = textLabel(fig.pan1,[0.05 0.20 0.4 0.05],...
    ['Min. Conn. = ' sprintf('%d',initConn)],fig.fS);
fig.edgeMinSlider = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.55 0.20 0.4 0.05],...
    'Style','slider',...
    'Min',0,...
    'Max',maxConn,...
    'Value',initConn,...
    'SliderStep',[5/maxConn 10/maxConn]);

% Activate slider
textLabel(fig.pan1,[0.05 0.15 0.4 0.05],'Scroll Wheel',fig.fS);
fig.activateSlider = uicontrol('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0.55 0.15 0.4 0.05],...
    'Style','popupmenu',...
    'String',{'None';'Transparency';'Connectivity'},...
    'Value',1,...
    'FontSize',fig.fS);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = textLabel(parent,pos,str,fS)
% Simple function to add a text label in the appropriate position

h = uicontrol('Parent',parent,....
    'Units','normalized',...
    'Position',pos,...
    'Style','text',...
    'String',str,...
    'FontSize',fS);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calculateNetworkXY(~,~,fig)
% Run the functions to draw the network

% Get the guidata
nw = guidata(fig.fig);

% Which method?
method = fig.method.String{fig.method.Value};

% Now perform
switch method
    
    case 'tSNE'
        
        % Get tSNE options
        dims = str2double(fig.tsneDims.String);
        perp = str2double(fig.tsnePerp.String);
        xy = tsne(nw.conn,[],2,dims,perp);
        
    case 'PCA'
        
        % Just run PCA on the matrix
        [~,xy,~] = pca(nw.conn,'NumComponents',2);
        
        
    case {'Force','Circle','Auto','Layered','Subspace'}
        tmpF = figure('Visible','off');
        hh = plot(nw.gr,'Layout',lower(method),'NodeLabel',[]);
        xy = [hh.XData' hh.YData'];
        close(tmpF);
        
    otherwise
        error('Otherwise');
end

% Just need to save the results to some kind of guidata, rather than
% drawing them in this function
nw.xy = xy;
guidata(fig.fig,nw);

% Now draw the results using the function
drawNetwork([],[],fig);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawNetwork(~,~,fig)
% Draw the network as selected...

tic

% Draw a waitbar to keep the user busy.
wb = waitbar(0.25,'Refreshing network');

% Get the guidata
nw = guidata(fig.fig);
if ~isfield(nw,'xy')
    return
end

% Reset the axes
axes(fig.ax); hold on;
f0 = get(fig.ax,'Children');
delete(f0);

% What would happen if we just used the default network plotting tool?
hh = plot(nw.gr,'Layout','force','NodeLabel',[]);

% Change the xy data to match...
hh.XData = nw.xy(:,1);
hh.YData = nw.xy(:,2);

% Let's just remove the node visibility by setting them to zero
hh.MarkerSize = 0.01;

% Do the same for the edges
hh.EdgeAlpha = 0;

% Save hh to guidata
nw.hh = hh;
guidata(fig.fig,nw);

% Change the edge / line properties
drawNetworkEdges([],[],fig);

% Change the node / scatter properties
drawNetworkNodes([],[],fig)

toc

delete(wb);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawNetworkEdges(~,~,fig)
% Change the properties of the lines connecting nodes

% Get the guidata
nw = guidata(fig.fig);

% Get values from the callbacks
minConn = fig.edgeMinSlider.Value;
maxWdth = str2double(fig.edgeMax.String);
tranVal = fig.transSlider.Value;

% Determine line widths
lw = nw.gr.Edges.Weight;

% Determine which line widths are below the threshold
belThr = lw < minConn;
lw = maxWdth * (lw / max(lw)) .^ 4;
lw(belThr) = 0.01;
nw.hh.LineWidth = lw;

% Set colours of edges to be white if low and black if high
cw = nw.gr.Edges.Weight;
cw = 90 * (cw / max(cw)) .^ 4;
cw = round(cw) + 1;
cw(belThr) = 1;
cmap = flipud(gray(101));
newc = cmap(cw,:);

% Update here
nw.hh.EdgeColor = newc;
nw.hh.EdgeAlpha = tranVal;

% Update guidata in case we need the updated settings
guidata(fig.fig,nw);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawNetworkNodes(~,~,fig)
% Change the nodes and scatter properties accordingly

% Get the guidata
nw = guidata(fig.fig);
if ~isfield(nw,'xy') || ~isfield(nw,'hh')
    return
end

% Reset the axes
axes(fig.ax); hold on;
f0 = get(fig.ax,'Children');
f0 = findobj(f0,'Type','scatter');
delete(f0);

% Which group to display as markers?
group = fig.marker.String{fig.marker.Value};

% Which group(s) to display the colours?
colList = fig.colour.String(fig.colour.Value);
if any(strcmp(colList,'NumAnno'))
    colList = {'NumAnno'};
    fig.colour.Value = 1;

end

% Determine colour groupings for the observations
cols = zeros(size(nw.xy,1),1);
for n = 1:numel(colList)
    tmp = nw.gr.Nodes.(colList{n});
    if isnumeric(tmp)
        cols = cols + tmp;
    end
end

% Now draw the locations...
mrk = {'o','d','s','p','h','^','v'};
grp = nw.gr.Nodes.(group);
[unq,~,ind] = unique(grp);
numG = numel(unq);
if isnumeric(unq)
    legLab = cell(numG,1);
    flag = true;
else
    flag = false;
end
lg = zeros(numG,1);
for n = 1:numG
    
    % Which ones?
    fx = ind == n;
    
    % Which marker?
    m = mod(n,numel(mrk));
    if m == 0
        m = numel(mrk);
    end
    
    % Draw the scatter
    ff = scatter(nw.xy(fx(1),1),nw.xy(fx(1),2),100,mrk{m},...
        'MarkerEdgeColor','k',...
        'Tag','scat');
    scatter(nw.xy(fx,1),nw.xy(fx,2),120,cols(fx,1),mrk{m},'filled',...
        'MarkerEdgeColor',[0.5 0.5 0.5],...
        'Tag','scat');
    
    lg(n,1) = ff(1);
    
    if flag
        legLab{n,1} = num2str(unq(n));
    else
        legLab{n,1} = unq{n};
    end
       
    
end

% Add the legend
legend(lg,legLab,...
    'FontSize',fig.fS,...
    'Location','South',...
    'Orientation','Horizontal');

% How many groups were selected?
if numel(colList) == 1
    cbLab = colList{1};
else
    cbLab = 'Multiple Group Presence/Absence';
end


% Colour bar formatting
cl = [min(cols) max(cols)];
if cl(1) == cl(2)
    cl(1) = cl(1) - 1;
end
caxis(cl);

if any(strcmp(colList,'NumAnno'))
    tickVals = 0:50:max(cols);
else
    tickVals = 0:1:max(cols);
end
set(fig.cb,'Ticks',tickVals);
ylabel(fig.cb,cbLab,...
    'FontWeight','bold',...
    'FontSize',fig.fS,...
    'Interpreter','none');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeEdgeTransparency(~,~,fig)
% Update the transparency on the graph...

% Get the guidata
nw = guidata(fig.fig);

% Get current value?
tval = fig.transSlider.Value;
fig.transValue.String = ['Trans = ' sprintf('%0.2f',tval)];

% Update the current values...
nw.hh.EdgeAlpha = tval;

% Update and move on
guidata(fig.fig,nw);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeEdgeMinConnectivity(~,~,fig)
% Update the minimum edge connectivity

% Get the guidata
nw = guidata(fig.fig);

% Get current value?
minConn = round(fig.edgeMinSlider.Value);
fig.edgeMinValue.String = ['Min. Conn. = ' sprintf('%d',minConn)];
tranVal = fig.transSlider.Value;

% Maximum line width
maxWdth = str2double(fig.edgeMax.String);

% Determine line widths
lw = nw.gr.Edges.Weight;

% Determine which line widths are below the threshold
belThr = lw < minConn;
lw = maxWdth * (lw / max(lw)) .^ 4;
lw(belThr) = 0.01;
nw.hh.LineWidth = lw;

% Set colours of edges to be white if low and black if high
cw = nw.gr.Edges.Weight;
cw = 90 * (cw / max(cw)) .^ 4;
cw = round(cw) + 1;
cw(belThr) = 1;
cmap = flipud(gray(101));
newc = cmap(cw,:);

% Update here
nw.hh.EdgeColor = newc;
nw.hh.EdgeAlpha = tranVal;

% Update and move on
guidata(fig.fig,nw);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sliderCallback(~,event,fig)
% Update the image in the gui

% Get the guidata and decide if we should continue
nw = guidata(fig.fig);
if ~isfield(nw,'xy')
    return
end

% We need to check which function the scroll wheel is to be controlling
ctrl = fig.activateSlider.String{fig.activateSlider.Value};
switch ctrl
    
    case 'None'
        % Just do nothing.
        return
        
    case 'Transparency'
        shiftVal = 0.05;
        handle = fig.transSlider;
        txtLab = 'Trans. = ';
        
    case 'Connectivity'
        shiftVal = 5;
        handle = fig.edgeMinSlider;
        txtLab = 'Min. Conn. = ';
end


% Which number to draw now?
val = handle.Value;

% What has happened?
switch event.EventName

    case 'WindowScrollWheel'            
        if event.VerticalScrollCount > 0
            
            val = val + shiftVal;
            
        elseif event.VerticalScrollCount < 0
            
            val = val - shiftVal;
            
        end

    case 'Action'  
        %disp('When does this action happen?');
        return
        
    otherwise
        disp('There is no otherwise');
        return
        
end
   
% Change the slider
if val < handle.Min
    val = handle.Min;
elseif val > handle.Max
    val = handle.Max;
end

% Now update...
%val = round(val);
handle.String = [txtLab sprintf('%d',val)];
handle.Value = val;

% Now update the edges between the nodes
switch ctrl
    case 'Transparency'
        changeEdgeTransparency([],[],fig);
    case 'Connectivity'
        changeEdgeMinConnectivity([],[],fig);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
