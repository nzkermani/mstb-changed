function stemPlotNice(mz,loadings,titText,ax)
% PLot the loadings

% Hold the axes
if nargin == 3
    ax = [];
end
if isempty(ax)
    figure;
else
    axes(ax);
end
hold on;

% Insert zeros...
%[x1,y1] = insertZeros(mz,loadings(:,1:2)',0.01);

% Apply an offset to make these visible
%y1(2,:) = y1(2,:) - (0.5 * max(y1(2,:))); 

%plot(mz,loadings(:,1:2));
stem(mz,loadings,'MarkerSize',0.01,'LineWidth',2);

box on;

title(titText,'FontSize',16);

set(gca,'YTickLabel',[],...
    'YTick',0,...
    'XTickMode','auto',...
    'XTickLabelMode','auto',...
    'LineWidth',5,...
    'TickLength',[0 0],...
    'YDir','normal');

grid off;

xlim([min(mz) max(mz)]);
ylim([min(loadings(:)) max(loadings(:))]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
