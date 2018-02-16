function statsTSNECalculate(src,~,fig,man)
%statsTSNECalculate - run the function to calculate the tsne, save it, and
%show it in the window

% Guidata
sts = guidata(fig.fig);

% Create a vector for intermediate visualisation
val = man.groups.Value;
str = man.groups.String;
grp = statsObservationLabels(sts.proc.meta,str,val);
[~,~,grpNum] = unique(grp);

% Let's check to see if we already have calculated the tsne values. In
% which case we jsut need to redraw it (rather than recalculate it).
if strcmp(src.Style,'listbox') && isfield(sts,'res')
    if isfield(sts.res,'tsne')        
        statsScatterPlot(fig.ax.scatter(1),sts.res.tsne.ss(:,1:2),grp,...
            false,false,sts.proc.meta);        
        return        
    end    
end

% Get the tsne options
opts.dims = str2double(man.dims.String);
opts.perp = str2double(man.perp.String);

% Set the axes as current to 'ensure' that the algorithm updates there...
szF = fig.fig.Position;
tmpFig = figure('Units','normalized',...
    'Position',szF);

% Run the function
[sts.res.tsne.ss] = tsne(sts.proc.sp,grpNum,2,opts.dims,opts.perp);
sts.res.tsne.opts = opts;
close(tmpFig);

% Save the results to guidata
guidata(fig.fig,sts);

% Draw the results in the correct axes
statsScatterPlot(fig.ax.scatter(1),sts.res.tsne.ss(:,1:2),grp,...
    false,false,sts.proc.meta);

% Ditch anything in the loadings axes
f0 = get(fig.ax.load(1),'Children');
delete(f0);
legend(fig.ax.load(1),'off');

% And remove old clusters etc.
f0 = get(fig.ax.conf(1),'Children');
delete(f0);
legend(fig.ax.conf(1),'off');



% And probably run other function parts as well...



end

