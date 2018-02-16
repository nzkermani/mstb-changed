function scatterPlotNice(scores,histID,cols,titText,ax)
% SCatter plot the MVA results

% Hold the axes
if nargin == 4
    ax = [];
end
if isempty(ax)
    figure;
else
    axes(ax);
end
hold on;

% Determine unique groups in the histIDs
[unq,~,ind] = unique(histID);
numG = numel(unq);
for n = 1:numG
    
    % Indices of points to plot
    fx = ind == n;
    
    symb = 'o';
    col = cols(n,:);
    
    
    % Scatter away
    scatter(scores(fx,1),scores(fx,2),120,col,symb,'filled',...
        'MarkerEdgeColor','k');
    
    % Determine the ellipse
    try
        [ell,~,~] = error_ellipse(scores(fx,:),95);
        plot(ell(:,1),ell(:,2),'Color',col,'LineWidth',2);
    catch
        disp('cannot determine ellipse');
    end
end

% Add a title...
title(titText,'FontSize',16);
    
box on;

set(gca,'XTickLabel',[],'YTickLabel',[],...
    'XTick',0,...
    'YTick',0,...
    'LineWidth',5,...
    'TickLength',[0 0]);
grid off;

% Draw lines along the origin
xlim('auto');
ylim('auto');
xl = xlim(gca);
yl = ylim(gca);
line([xl(1) xl(2)],[0 0],'LineStyle',':','Color','k');
line([0 0],[yl(1) yl(2)],'LineStyle',':','Color','k');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
