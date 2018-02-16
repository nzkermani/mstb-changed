function [ output_args ] = imzmlRecalPlot(devs,files,timeStamp)
%imzmlRecalPlot - create a nice fancy plot

figure('Units','pixels','Position',[286 890 1274 448]); hold on;

med = nanmedian(devs,1);

[~,ord] = sort(med);

boxplot(devs(:,ord));

% Set y limits
yl = ylim;
ylim([-max(abs(yl)) max(abs(yl))]);

% Determine x limits
xl = size(devs,2);
xlim([0 xl+1])

% Draw a line
line([0 xl+1 NaN 0 xl+1 NaN 0 xl+1],[3 3 NaN -3 -3 NaN 0 0]);

% Remove the ticks...
set(gca,'XTick',[]);
set(gca,'FontSize',14);
ylabel('ppm deviation','FontSize',16,'FontWeight','bold');


% Make another figure using the timeStamp
figure('Units','pixels','Position',[286 890 1274 448]); hold on;

scatter(timeStamp,med,80,'r','o','filled','MarkerEdgeColor','k');

% Mark those which are outside of the limits
ol = find(abs(med) > 3);
text(timeStamp(ol)+10,med(ol),files(ol,1),...
    'FontSize',12,...
    'Rotation',20);

xl = xlim;

line([xl NaN xl NaN xl],[3 3 NaN -3 -3 NaN 0 0]);

ylim([-max(abs(yl)) max(abs(yl))]);

xlim(xl);

box on;

datetick('x','yy-mm');

set(gca,'FontSize',14);
xlabel('Date','FontSize',16,'FontWeight','bold');
ylabel('Median ppm deviation','FontSize',16,'FontWeight','bold');

end

