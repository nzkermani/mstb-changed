function defGraphLabs(xLab,yLab,cBar)
%defGraphLabs - apply a consistent set of labels and fontsizes to the
%current figure

% Define defaults
opts.axes.size  = 14;

opts.label.size = 16;
opts.label.wght = 'bold';
opts.label.angl = 'normal';

opts.cbar.cmap  = 'jet';

% Auto-format m/z into italics
if strcmp(xLab,'m/z')
    opts.label.angl = 'italic';
end

% X axis labels
xlabel(xLab,...
    'FontSize',opts.label.size,...
    'FontWeight',opts.label.wght,...
    'FontAngle',opts.label.angl);

% Y axis labels
ylabel(yLab,...
    'FontSize',opts.label.size,...
    'FontWeight',opts.label.wght);

% Optional colour bar addition...
if nargin == 3
    c = colorbar;
    colormap(opts.cbar.cmap);
    ylabel(c,cBar,...
        'FontSize',opts.label.size,...
        'FontWeight',opts.label.wght);
end

% Axes properties
box on;
set(gca,'FontSize',opts.axes.size);

set(gcf,'Position',[400 400 1200 800]);



end

