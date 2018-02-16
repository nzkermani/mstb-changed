function statsLDADraw(~,~,fig,window)
%statsLDADraw - draw the plots from LDA

% Guidata
sts = guidata(fig.fig);

% Determine the components...
try
    cs = [get(window.comp1,'Value') get(window.comp2,'Value')];
catch
    cs = [1 2];
end

% Plot the ellipses?
try
    doEllipse = get(window.ellipse,'Value');
catch
    doEllipse = false;
end

% Work out what colour group is selected for colouring in. This is to allow
% file biases to be visualised, rather than anything else.  ALthough it
% does have the potential for abuse...
if isfield(window,'colgroup')
    str = window.colgroup.String;
    val = window.colgroup.Value;
    metaGrp = statsObservationLabels(sts.proc.meta,str,val);    
else
    metaGrp = sts.res.mmc.grp;
end
    

% Scores
statsScatterPlot(fig.ax.scatter(1),sts.res.mmc.ss(:,cs),metaGrp,doEllipse,false,sts.proc.meta);

% Loadings

% Loadings
if strcmp(sts.datatype,'ms')
    
    statsLoadingsPlot(fig.ax.load(1),...
        sts.proc.var.mz,...
        sts.res.mmc.ll(:,cs));
    
elseif strcmp(sts.datatype,'lcms')    
    statsLoadingsPlotScatter(fig.ax.load(1),...
        sts.proc.var.mz,...
        sts.proc.var.rt,...
        sts.res.mmc.ll(:,cs));
    
    % What are the axes limits of the plot next door?
    %axlx = xlim(fig.ax.spec);
    %axly = ylim(fig.ax.spec);
    %xlim(fig.ax.load(1),axlx);
    %ylim(fig.ax.load(1),axly);
    
end


% The extra plot
try
    statsLDAExtra([],[],fig,window);
catch
    parent = fig.ax.conf(1);
    f0 = get(parent,'Children');
    delete(f0);
end

% Link axes?
if strcmp(sts.datatype,'ms')
    linkaxes([fig.ax.load(1) fig.ax.spec(1)],'x');
else
    linkaxes([fig.ax.load(1) fig.ax.spec(1)],'xy');
end

end

