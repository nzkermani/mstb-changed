function statsPlotSpectra(~,~,fig,type)
%statsPlotSpectra - draw the spectra in the plot...

sts = guidata(fig.fig);
if isempty(sts)
    return
end

% Simply define the parent axes in which this going to be drawn
parent = fig.ax.spec(1);
f0 = get(parent,'Children');
delete(f0);

axes(parent);
hold on;

% What kind of data to plot? The raw or the processed?
switch sts.datatype
    
    case 'ms'
        plotMS(sts,type);
        
    case 'lcms'
        plotLCMS(sts,type);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotMS(sts,type)

switch type
    case 'raw'
        mz = sts.raw.var.mz;
        sp = sts.raw.sp;
        
    case 'proc'
        mz = sts.proc.var.mz;
        sp = sts.proc.sp;
end

% Insert zeros...
[x1,y1] = insertZeros(mz,sp,0.01);

% Simple plot of all together
plot(x1,y1);
set(gca,'YTickLabel',[],...
    'YTick',0,...
    'XTickMode','auto',...
    'XTickLabelMode','auto',...
    'LineWidth',5,...
    'TickLength',[0 0],...
    'YDir','normal',...
    'FontSize',14,...
    'FontWeight','bold');

grid off;

xlim([min(x1) max(x1)]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotLCMS(sts,type)

switch type
    case 'raw'
        mz = sts.raw.var.mz;
        rt = sts.raw.var.rt;
        sp = sts.raw.sp;        
        
    case 'proc'
        mz = sts.proc.var.mz;
        rt = sts.proc.var.rt;
        sp = sts.proc.sp;
end

% Scatter plot it all - using the mean/median intensity for colour...
sp = nanmean(sp,1);

scatter(mz,rt,20,sp,'o','filled','MarkerEdgeColor','k');

set(gca,...
    'YTickMode','auto',...
    'YTickLabelMode','auto',...
    'XTickMode','auto',...
    'XTickLabelMode','auto',...
    'LineWidth',5,...
    'TickLength',[0 0],...
    'YDir','normal',...
    'FontSize',12,...
    'FontWeight','bold');

% set(gca,'YTickLabel',[],...
%     'YTick',0,...
%     'XTickMode','auto',...
%     'XTickLabelMode','auto',...
%     'LineWidth',5,...
%     'TickLength',[0 0],...
%     'YDir','normal',...
%     'FontSize',12,...
%     'FontWeight','bold');

grid off;

xlim([min(mz) max(mz)]);
ylim([min(rt) max(rt)]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
