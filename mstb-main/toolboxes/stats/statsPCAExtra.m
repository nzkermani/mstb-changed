function statsPCAExtra(parent,type,mz,res,cs,grp)
% statsPCAExtra - plot some things in the extra axes

% Delete existing objects
f0 = get(parent,'Children');
delete(f0);
legend(parent,'off');

% Hold the axes
axes(parent);
hold on;

% Switch through the options...
switch type
    case 'Loadings'
        plotLoads(mz,res.ll(:,cs));
        
    case 'Eigenvalues'
        plotEigs(res.ee);
        
    case 'Boxplots'
        % Plot something
        disp('Plot something');
        
end

set(gca,'XTickLabel',[],...
    'YTickLabel',[],...
    'XTick',0,...
    'YTick',0,...
    'LineWidth',5,...
    'TickLength',[0 0],...
    'YDir','normal');

grid off;

return

% Draw lines along the origin
xlim('auto');
ylim('auto');
xl = xlim(gca);
yl = ylim(gca);
line([xl(1) xl(2)],[0 0],'LineStyle',':','Color','k');
line([0 0],[yl(1) yl(2)],'LineStyle',':','Color','k');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotEigs(ee)

maxE = min([10 numel(ee)]);

stem(1:maxE,ee(1:maxE),...
    'Color','k',...
    'MarkerFaceColor','k',...
    'MarkerEdgeColor','k',...
    'LineWidth',4,...
    'MarkerSize',10);

txt = int2str(ee(1:maxE));

text((1:maxE)+0.2,ee(1:maxE),txt,...
    'Rotation',00,...
    'FontSize',18);

xlim([0 maxE+1]);
ylim([-1 ee(1)+5]);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotLoads(mz,ll)

% Find the 5 loadings the furthest from zero...
dist = sqrt(sum(ll .^ 2,2));
[~,idx] = sort(dist,'descend');
idx = idx(1:5);

scatter(ll(:,1),ll(:,2),80,mz,'o','filled',...
    'MarkerEdgeColor','k');

labs = cell(size(idx));
for n = 1:size(labs,1)
    labs{n,1} = ['  ' sprintf('%0.2f',mz(idx(n)))];
end

text(ll(idx(1:5),1),ll(idx(1:5),2),labs,...
    'FontSize',14)

xl = [min(ll(:,1)) max(ll(:,1))] * 1.2;
yl = [min(ll(:,2)) max(ll(:,2))] * 1.2;

xlim(xl);
ylim(yl);

line([xl(1) xl(2)],[0 0],'LineStyle','-','Color',[16 135 232]/256);
line([0 0],[yl(1) yl(2)],'LineStyle','-','Color',[230 115 142]/256);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
