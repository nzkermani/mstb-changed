function statsLoadingsPlotScatter(parent,mz,rt,loadings)
% PLot the loadings

f0 = get(parent,'Children');
delete(f0);

axes(parent);
hold on;

%scatter(loadings(:,1),loadings(:,2),80,mz,'o','filled');

% Determine size based on distance from the origin - bigger is better
dx = loadings .^ 2;
dx = sqrt(sum(dx,2));
dx = dx .^ 4;

% Scale accordingly... biggest of say 120
dx = 150 * dx / max(dx);
fx = dx > 1;
disp(['Displaying ' int2str(sum(fx)) '/' int2str(numel(fx)) ' variables']);

scatter(mz(fx),rt(fx),dx(fx),dx(fx),'o','filled');

set(gca,'XTickLabelMode','auto',...
    'XTickMode','auto',...
    'YTickLabelMode','auto',...
    'YTickLabel','auto',...
    'LineWidth',5,...
    'TickLength',[0 0],...
    'FontSize',12,...
    'FontWeight','bold',...
    'YDir','normal');

grid off;

% xlim([min(loadings(:,1)) max(loadings(:,1))]);
% ylim([min(loadings(:,2)) max(loadings(:,2))]);

xlim([min(mz) max(mz)]);
ylim([min(rt) max(rt)]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
