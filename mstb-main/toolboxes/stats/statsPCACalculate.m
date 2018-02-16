function statsPCACalculate(src,event,fig,window)
%statsPCACalculate - get the bits required to calculate or just recolour
%the PCA scores plot...

wb = waitbar(0.5,'PCA');

% Guidata
sts = guidata(fig.fig);

% First thing is to see if we have already calculated it. If not, then we
% need to do so immediately...
if ~isfield(sts.res,'pca')
    flag = true;
elseif isempty(sts.res.pca)
    flag = true;
else
    flag = false;
end

% Calculate
if flag    
    [sts.res.pca.ll,sts.res.pca.ss,ee] = pca(sts.proc.sp,'NumComponents',10);
    sts.res.pca.ee = 100 * ee / sum(ee);
    
    maxComp = size(sts.res.pca.ss,2);
    set(window.comp1,'Value',1,'String',int2str([1:maxComp]'));
    set(window.comp2,'Value',2,'String',int2str([1:maxComp]'));
    
end

% Get the selected value of metadata from the table...
val = window.groups.Value;
str = window.groups.String;
grp = statsObservationLabels(sts.proc.meta,str,val);

% Determine the components...
cs = [get(window.comp1,'Value') get(window.comp2,'Value')];

% Plot the ellipses?
doEllipse = get(window.ellipse,'Value');

% Plot centroids?
doCentroids = get(window.centroid,'Value');

% Scores
statsScatterPlot(fig.ax.scatter(1),sts.res.pca.ss(:,cs),grp,...
    doEllipse,doCentroids,sts.proc.meta);

% Loadings
if strcmp(sts.datatype,'ms')    
    statsLoadingsPlot(fig.ax.load(1),...
        sts.proc.var.mz,...
        sts.res.pca.ll(:,cs));

elseif strcmp(sts.datatype,'lcms')    
    statsLoadingsPlotScatter(fig.ax.load(1),...
        sts.proc.var.mz,...
        sts.proc.var.rt,...
        sts.res.pca.ll(:,cs));
    
    % What are the axes limits of the plot next door?
    axlx = xlim(fig.ax.spec);
    axly = ylim(fig.ax.spec);
    xlim(fig.ax.load(1),axlx);
    ylim(fig.ax.load(1),axly);
    
end

% The extra plot
val = window.extraPlot.Value;
str = window.extraPlot.String;
statsPCAExtra(fig.ax.conf(1),str{val},sts.proc.var.mz,sts.res.pca,cs,grp);


% Link axes?
if fig.ax.spec.Position(3) > 0.95
    linkaxes([fig.ax.load(1) fig.ax.spec(1)],'x');
else
    linkaxes([fig.ax.load(1) fig.ax.spec(1)],'xy');
end

% Update GUI
guidata(fig.fig,sts);

% Delete the waitbar, although it appears to go before the plot has been
% properly updated
delete(wb);

end

