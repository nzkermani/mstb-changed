function statsUnivariateVolcano(~,~,fig,window)
% Stats UV volcano

% Get the guidata
sts = guidata(fig.fig);
if ~isfield(sts,'res') || ~isfield(sts.res,'uv')
    return
end

% Where are we drawing the things
parent = fig.ax.scatter;
parent2 = fig.ax.conf;

% What are the colours?
unqG = unique(sts.res.uv.grp);
numG = numel(unqG);
cols = parula(numG);

% Determine useful thresholds - eventually these will be defined by the
% user...
pqThresh = -log10(str2double(get(window.pqThresh,'String')));
fcThresh = str2double(get(window.fcThresh,'String'));

% Plot p or q values?
pqChoice = get(window.plotPQ,'Value');
if pqChoice == 1
    pqString = 'p';
else
    pqString = 'q';
end

% Log the pq values and use these throughout.  Need to remove the Inf
% values, or set them to be above the others by a slight way
pq = sts.res.uv.pq(:,pqChoice);
fc = sts.res.uv.fc;
lowPQ = min(pq(pq > 0));
idxLo = pq == 0;
pq(isnan(pq)) = 1;
pq = -log10(pq);
pq(idxLo) = -log10(lowPQ) * 1.2;

% Do below but for the plot above
f0 = get(parent2,'Children');
title(parent2,'');
delete(f0);

% Delete existing objects
f0 = get(parent,'Children');
delete(f0);
legend(parent,'off')

% Hold the axes
axes(parent);
hold on;

% Highlight only significant ones
fx = abs(fc) >= fcThresh & pq >= pqThresh;
if sum(fx) == 0
    fx = false(size(pq));
    areSig = false;
else
    areSig = true;
end

% Determine X and Y axes limits
absFC = max(abs(fc)) * 1.1;
if absFC <= fcThresh
    xl = [-fcThresh-0.2 fcThresh+0.2];
else
    xl = [-absFC absFC];
end
yl = [-0.1 max(pq) * 1.05];

% We need to draw the gridlines before the other points, as we don't want
% the gridlines to be on top of the points
line([xl(1) xl(2)],[pqThresh pqThresh],...
    'LineStyle','--',...
    'Color','k');

% Add in another for a higher value...
line([xl(1) xl(2)],[floor(max(pq)) floor(max(pq))],...
    'LineStyle',':',...
    'Color','k');

% Vertical lines for FC significance
line([-fcThresh -fcThresh],[yl(1) yl(2)],...
    'LineStyle','--',...
    'Color','k');
line([fcThresh fcThresh],[yl(1) yl(2)],...
    'LineStyle','--',...
    'Color','k');

% Define colours
colsM = [0.2627 0.5765 0.7647; 0.8392 0.3765 0.3020];
% Scatter away - bad variables
scatter(fc(~fx),pq(~fx),60,colsM(1,:),'o','filled',...
    'MarkerEdgeColor',[0.7 0.7 0.7]);
    
% These are the interesting ones
switch sts.datatype
    case 'ms'
        try
            xy = [fc(fx) pq(fx) sts.proc.var.mz(fx)];
        catch
            xy = [fc(fx) pq(fx) sts.proc.var.mz(fx)'];
            sts.proc.var.mz = sts.proc.var.mz';
        end
    case 'lcms'
        xy = [fc(fx) pq(fx) sts.proc.var.mz(fx) sts.proc.var.rt(fx)];
end

% Only try plotting significant features if there are any...
if areSig
    scatter(fc(fx),pq(fx),120,colsM(2,:),'o','filled',...
        'MarkerEdgeColor',[0.7 0.7 0.7],...
        'HitTest','on',...
        'ButtonDownFcn',{@statsUnivariateClick,...
        [parent parent2],...
        fc(fx),pq(fx),...
        sts.proc.sp(:,fx),...
        sts.res.uv.grp,...
        xy,...
        cols,...
        pqString});
end

% Add text labels so that we can make sure that the boxplots are shown
% correctly
%text(fc(fx),pq(fx),num2str(mz(fx)'));
    
box on;

% Try to determine the Y axis limit
if sum(fx) == 0
    yTickLim = pqThresh;
    yTickLab = {['q = ' num2str(10^-pqThresh)]};
else
    yTickLim = [pqThresh floor(max(pq))];
    yTickLab = {['q = ' num2str(10^-pqThresh)],...
        ['q = ' num2str(10^-yTickLim(2))]};
end

set(gca,...
    'XTickLabel',[-fcThresh fcThresh],...
    'YTickLabel',yTickLab,...
    'XTick',[-fcThresh fcThresh],...
    'YTick',yTickLim,...
    'LineWidth',5,...
    'TickLength',[0 0],...
    'YDir','normal');
grid off;

xlim(xl);
ylim(yl);


% Display the significant ones
if areSig
    ff = array2table([sts.proc.var.mz(fx) pq(fx) fc(fx)],...
        'VariableNames',{'m_z','neg_log10_q','log2_FC'});
    disp(ff);
else
    disp('No significant variables');
    disp('Try lowering the thresholds');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

