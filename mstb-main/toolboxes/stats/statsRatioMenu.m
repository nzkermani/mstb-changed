function statsRatioMenu(~,~,fig)
%statsRatioMenu

% Get the guidata
sts = guidata(fig.fig);
if isempty(sts)
    return
end

% Determine the correlation matrix, unless already calculated
if ~isfield(sts.res,'ratio')
    
    % Display waitbar so that we know something is happening
    wb = waitbar(0.5,'Calculating correlations');
    
    % Calculate correlations!
    [a,b] = corr(sts.proc.sp);
    
    % Set bottom part to NaN
    low = tril(ones(size(a)),0);
    a(low == 1) = NaN;
    b(low == 1) = NaN;

    % Set poor pvals to NaN
    mask = b > 0.01;
    a(mask) = NaN;
    
    % Save
    sts.res.ratio.rval = a;
    sts.res.ratio.pval = b;

    % Set poor pval to NaN
    guidata(fig.fig,sts);
    
    % Delete the waitbar like a good boy
    delete(wb);
end

% Draw the window to help control the colours / groups / etc
[window] = manipulateWindow(fig,sts);

% Plot the correlations...
plotCorrelations([],[],fig,window);

set(window.threshold,'Callback',{@plotCorrelations,fig,window});
set(window.var1,'Callback',{@variableBoxClick,fig,window});
set(window.var2,'Callback',{@variableBoxClick,fig,window});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [man] = manipulateWindow(fig,sts)
% Window with the little options for manipulating the figure...

% This is where we draw everything
parent = fig.pan2;

fS = 14;

% Heading
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.95 1 0.05],...
    'Style','text',...
    'String','Ratio Analyses',...
    'FontSize',24,...
    'BackgroundColor',[1 1 1]);

%%%%

% Correlation threshold
x = 0.9;
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 x 0.5 0.05],...
    'Style','text',...
    'String','Threshold (negative)',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
man.threshold = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.5 x+0.025 0.475 0.02],...
    'Style','edit',...
    'String','0.8',...
    'FontSize',fS);%'Min',1,...'Max',3);


% Heading
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 0.7 1 0.05],...
    'Style','text',...
    'String','Visualisation',...
    'FontSize',18,...
    'BackgroundColor',[1 1 1]);

% Which meta group for boxplot visualisation
x = 0.65;
fn = fieldnames(sts.proc.meta);
if any(strcmp(fn,'histID'))
    fnval = find(strcmp(fn,'histID'));
else
    fnval = 1;
end    
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 x 0.33 0.05],...
    'Style','text',...
    'String','Group',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
man.groupSelect = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.34 x-0.15 0.64 0.2],...
    'Style','popupmenu',...
    'String',fn,...
    'Value',fnval,...
    'FontSize',fS);

% Variables worth plotting
x = 0.6;
varN = sprintf('%0.3f\n',sts.proc.var.mz);
uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0 x 0.33 0.05],...
    'Style','text',...
    'String','Variables',...
    'FontSize',fS,...
    'BackgroundColor',[1 1 1]);
man.var1 = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.33 x-0.15 0.32 0.2],...
    'Style','listbox',...
    'String',varN,...
    'Value',1,...
    'FontSize',fS,...
    'BackgroundColor',[0.9 0.9 0.9],...
    'ForegroundColor','k');%'Min',1,...'Max',3);
man.var2 = uicontrol('Parent',parent,...
    'Units','normalized',...
    'Position',[0.66 x-0.15 0.32 0.2],...
    'Style','listbox',...
    'String',varN,...
    'Value',1,...
    'FontSize',fS,...
    'BackgroundColor',[0.9 0.9 0.9],...
    'ForegroundColor','k');%'Min',1,...'Max',3);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotCorrelations(src,event,fig,man)
% Scatter plot significant correlations

% Get the threshold - convert to negative
thresh = str2double(man.threshold.String) * -1;

% Get corr matrix
sts = guidata(fig.fig);
rval = sts.res.ratio.rval;
%pval = sts.res.ratio.pval;

% Extract the good ones, for scatter plot analysis
% Find all values of rval which are true...
[fx,fy] = find(rval < thresh);

% Warning for too many points to handle
if numel(fx) > 10000
    wd = warndlg('Too many points to display');
    pause(1);
    delete(wd);
    return
end

% Prep the data
data = zeros(numel(fx),3);
for n = 1:numel(fx)    
    data(n,1:3) = [sts.proc.var.mz(fx(n)) sts.proc.var.mz(fy(n)) rval(fx(n),fy(n))];
end

% Axes prep.
f0 = get(fig.ax.scatter,'Children');
delete(f0);
axes(fig.ax.scatter);
hold on;
legend(fig.ax.scatter,'off')

% Set the axes limits...
xl = [min(sts.proc.var.mz)-10 max(sts.proc.var.mz)+10];
yl = xl;
xlim(fig.ax.scatter,xl);
ylim(fig.ax.scatter,yl);

% Draw some lines and useful information...
line(xl,yl,'Color','r','LineWidth',2,'LineStyle','-.');

% Custom gridlines, every 100 m/z
xg = 0:100:xl(2);
xg = [zeros(size(xg)); xg; NaN(size(xg))];
yg = [xg(2,:); xg(2,:); xg(3,:)];

yg2 = [xg(2,:); xg(2,:); xg(3,:)];
yg2(2,:) = max(xg(:,end));

% Plot them...
lcol = [0.7 0.7 0.7];
line(xg,yg,'Color',lcol,'LineWidth',0.5,'LineStyle',':');
line(yg,yg2,'Color',lcol,'LineWidth',0.5,'LineStyle',':');


% Scatter plot with appropriate click function...
scatter(fig.ax.scatter,...
    data(:,1),...
    data(:,2),...
    80,....
    data(:,3),...
    'o','filled',...
    'ButtonDownFcn',{@ratioClick,fig,man,data});

% Text labels at approp. coords. - try to ignore the first and last ones as
% we don't want labels spilling off the axes
for n = 2:size(yg,2)-0    
    text(yg(1,n)+10,...
        yg(2,n)-10,...
        int2str(yg(2,n)),...
        'FontSize',12,...
        'Clipping','on',...
        'Rotation',-45);    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ratioClick(~,event,fig,man,data)
% Clicked callback function for the main correlation image

% Get the guidata
sts = guidata(fig.fig);

% Get the click points
cx = event.IntersectionPoint(1);
cy = event.IntersectionPoint(2);

% Determine closest in 'data'
df = [data(:,1)-cx data(:,2)-cy];
df = df .^ 2;
[~,idx] = min(sum(df,2));

% Determine indices in mz
[~,ix] = min(abs(sts.proc.var.mz - data(idx,1)));
[~,iy] = min(abs(sts.proc.var.mz - data(idx,2)));

% Determine the ratio of these two variables...
rto = sts.proc.sp(:,ix) ./ sts.proc.sp(:,iy);
rto = log2(rto);
rto(isinf(rto)) = NaN;
%fx = ~isnan(rto);

% Determine which meta group
str = man.groupSelect.String;
val = man.groupSelect.Value;
metaGrp = statsObservationLabels(sts.proc.meta,str,val);

% Set the values of the boxes in man.val1,val2 to reflect the current
% selection
man.var1.Value = ix;
man.var2.Value = iy;

% Run the actual plotting function now that we have worked everything
% out...
doActualPlotting(fig.ax.conf,rto,metaGrp,data(idx,:));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function variableBoxClick(~,~,fig,man)
% What happens when we click the box in the side menu

% Guidata
sts = guidata(fig.fig);

% Determine the two variables...
v1 = man.var1.Value;
v2 = man.var2.Value;

% COnvert to mz to check
m1 = sts.proc.var.mz(v1);
m2 = sts.proc.var.mz(v2);

% Determine which meta group
str = man.groupSelect.String;
val = man.groupSelect.Value;
metaGrp = statsObservationLabels(sts.proc.meta,str,val);

% Correlate the variables...
crr = corr(sts.proc.sp(:,v1),sts.proc.sp(:,v2));

% Ratio the variables
rto = sts.proc.sp(:,v1) ./ sts.proc.sp(:,v2);
rto = log2(rto);
rto(isinf(rto)) = NaN;

% Now plot the results
doActualPlotting(fig.ax.conf,rto,metaGrp,[m1 m2 crr]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doActualPlotting(parent,rto,metaGrp,data)

% Plot the boxplot in the neighbouring figure...
f0 = get(parent,'Children');
delete(f0);
axes(parent);
hold on;
legend(parent,'off');

% If there is only one kind of value, then nothing to do here
if numel(unique(rto)) == 1
    xlim(parent,[0 1]);
    ylim(parent,[0 1]);
    set(parent','XTick',[],'YTick',[]);
    text(0.2,0.5,'Ratios all zero','FontSize',30);

    return
end

% Determine unique bits about meta
[unq,~,~] = unique(metaGrp);
numG = numel(unq);
cols = parula(numG);

% Use my default boxplot function
jsmBoxPlot(rto,metaGrp,...
    'Orientation',parent,...
    'Colours',cols,...
    'Legend',false,...
    'Labels',unq,...
    'Order',1:numel(unq));

% Change some of the axes properties
set(parent,...
    'YTickLabel',[],...
    'FontSize',10,...
    'FontWeight','normal',...
    'FontName','Helvetica',...
    'YDir','normal',...
    'YTickMode','auto',...
    'YTickLabelMode','auto');

ylabel(parent,'');

% Add some text in the bottom to tell us what the variables are...
yl = ylim(parent);

span = (yl(2) - yl(1)) * 1.1;
textY = yl(1) + ((yl(2) - yl(1)) * 1.055);
yl(2) = yl(1) + span;
ylim(parent,yl);

txt = ['Log_2 ratios of m/z = ' sprintf('%0.2f',data(1)) ' / ' ...
    sprintf('%0.2f',data(2)) ...
    ', r = ' sprintf('%0.2f',data(3))];
textX = 0.6;
text(textX,textY,txt,'FontSize',14);

line([0 numel(unq)+0.5],[0 0],...
    'LineWidth',2,...
    'LineStyle',':',...
    'Color',[0.7 0.7 0.7]);



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
