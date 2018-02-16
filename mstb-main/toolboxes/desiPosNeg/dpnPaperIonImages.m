function dpnPaperIonImages(extr,totsum,crd)
% Make image plots for showing poor ions in the dpnCorrelation plots

% How many ions?
numI = size(extr,3);

% Make a figure
figure('Units','normalized','Position',[0 0.25 1 0.5]);

% Axes handles
ax = zeros(1,numI+1);

% Loop through each ion
for i = 1:numI

    ax(i) = subplot(1,numI+1,i); hold on;
    scatter([0.5 crd(2)],[crd(1) 0.5],80,'k','d','filled');
    set(gca,'YDir','reverse');
    imagesc(extr(:,:,i));
    axis square;
    axis off;
    set(gca,'FontSize',16);
    hold on;
    
end

% Draw the total ion image at the end of the run
i = numI + 1;
ax(i) = subplot(1,numI+1,numI+1); hold on;
scatter([0.5 crd(2)],[crd(1) 0.5],80,'k','d','filled');
set(gca,'YDir','reverse');
imagesc(totsum);
axis square;
axis off;
set(gca,'FontSize',16);
hold on;

%  Link the axes
linkaxes(ax,'xy');

end