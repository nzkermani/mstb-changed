function statsBoxplotsDraw(~,~,fig,man)
%statsBoxplotsDraw - draw boxplots based on the selection

% Get the guidata
sts = guidata(fig.fig);

% Get the value of the clicked items...
str = man.list.String;
val = man.list.Value;
mzrt = str{val};

% And the groupings
str = man.groups.String;
val = man.groups.Value;
[grp] = statsObservationLabels(sts.proc.meta,str,val);

% MZRTRTMZ
format = man.mzrtchoice.Value;

% What about where to plot it?
newFig = man.external.Value;

% What about ROC instead?
plotROC = man.plotROC.Value;

% Separate mz | rt into separate parts, dependent on which comes first
fx = strfind(mzrt,' | ');
switch format
    case 1
        mz = str2double(mzrt(1:fx(1)-1));
        rt = str2double(mzrt(fx(1)+3:end));
    case 2
        rt = str2double(mzrt(1:fx(1)-1));
        mz = str2double(mzrt(fx(1)+3:end));
end

% So now we need to find the best matching of the variables...
dmz = (sts.proc.var.mz - mz) .^ 2;
%drt = (sts.proc.var.rt - rt) .^ 2;
drt = (sts.proc.var.mz - mz) .^ 2;

ddd = dmz + drt;
[~,b] = min(ddd);

% Unique group ordering
[unq,~,ind] = unique(grp);

% Plot a ROC instead?
if numel(unq) == 2 && plotROC
    
    % Delete existing
    f0 = findobj('Tag','ROCCURVE');
    close(f0);
    
    meanVal = [mean(sts.proc.sp(ind == 1,b)) mean(sts.proc.sp(ind == 2,b))];
    [~,rocIdx] = max(meanVal);
    
    [roc.a,roc.b,roc.c,roc.d,roc.e] = perfcurve(grp,sts.proc.sp(:,b),unq{rocIdx});
    
    % Find optimum values...
    opt1 = roc.a == roc.e(1);
    opt2 = roc.b == roc.e(2);
    optX = opt1 & opt2;
    optThr = roc.c(optX);
    
    figure('Units','pixels','Position',[600 612 974 420],'Tag','ROCCURVE');
    ax1 = subplot(1,2,1); hold on;
    randY = rand(numel(grp),1);
    r1 = ind == 1;
    r2 = ind == 2;
    rand1 = rand(sum(r1),1);
    rand2 = rand(sum(r2),1) + 1;
    lg1 = scatter(sts.proc.sp(r1,b),rand1,'red','o','filled');
    lg2 = scatter(sts.proc.sp(r2,b),rand2,'blue','o','filled');
    plot([optThr optThr],[0 2.5],'LineWidth',2,'Color','k');
    legend([lg1(1) lg2(1)],unq,'Location','NorthWest')
    box on;
    axis square;
    title(['m/z = ' sprintf('%0.3f',mz)]);
    xlabel('Intensity');
    ylabel('Random Index');
    set(gca,'FontSize',14,'YTick',[]);
    
    ax2 = subplot(1,2,2); hold on;
    plot(roc.a,roc.b,'LineWidth',2,'Color','k');
    scatter(roc.a(opt1&opt2),roc.b(opt1&opt2),80,'r','o','filled');
    box on;
    axis square;
    xlabel('FDR');
    ylabel('TPR');
    title(['ROC Curve; AUC = ' sprintf('%0.2f',roc.d)]);
    set(gca,'FontSize',14);
    
    return
else
    disp('ROC Not possible');
    man.plotROC.Value = 0;
end
    

% Determine a p/q value for
[pval,~,~] = anova1(sts.proc.sp(:,b),grp,'off');
[pval2,~,~] = kruskalwallis(sts.proc.sp(:,b),grp,'off');

disp(['ANOVA p = ' sprintf('%0.3E',pval)]);
disp(['KRUSKAL p = ' sprintf('%0.3E',pval2)]);

% Define the colours...
cols = parula(numel(unq));

% Reset the secondary axes
switch newFig
    case 0
        parent = fig.ax.conf;
        axes(parent);
        f0 = get(parent,'Children');
        delete(f0);
    case 1
        %nf = figure('Units','normalized','Position',[0.25 0.25 0.5 0.5]);
        
        jsmBoxPlot(sts.proc.sp(:,b),ind,...
            'Orientation','vertical',...
            'Colours',cols,...
            'Legend',false,...
            'Labels',unq,...
            'Order',1:numel(unq));        
        return
        
end
        
% Draw the boxplot in the figure axes
bpmethod = 'jsm';
switch bpmethod
    
    case 'jsm'

        % Use my default boxplot function
        jsmBoxPlot(sts.proc.sp(:,b),ind,...
            'Orientation',parent,...
            'Colours',cols,...
            'Legend',false,...
            'Labels',unq,...
            'Order',1:numel(unq));
       
        % Change some of the axes properties
        set(fig.ax.conf,...
            'YTickLabel',[],...
            'FontSize',10,...
            'FontWeight','normal',...
            'FontName','Helvetica',...
            'YDir','normal',...
            'YTickMode','auto',...
            'YTickLabelMode','auto');
        
        ylabel(parent,'');
        %title(fig(2),'');
        
        % Determine axes limits
        yl = [min(sts.proc.sp(:,b)) max(sts.proc.sp(:,b))];
        rnge = yl(2) - yl(1);        
        if rnge == 0
            rnge = 1;
        end
        yl(2) = yl(2) + (rnge*0.1);
        ylim(fig.ax.conf,yl);
        
        txt = ['m/z = ' sprintf('%0.3f',mz) ', p = ' sprintf('%0.4E',pval)];
        text(0.6,yl(2) - (rnge * 0.05),txt,...
            'FontSize',18,'FontWeight','bold',...
            'HorizontalAlignment','left');

    case 'trad'

        % Now do the boxplot of this variable...
        boxplot(parent,...
            sts.proc.sp(:,dd),...
            ind,... % groups
            'Labels',unq,...
            'Jitter',0.8);

end



end

