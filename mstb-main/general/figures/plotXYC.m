function [fig,ax] = plotXYC(x,y,col,cbLab)
%plotXYC - colour the y data using the colours in col
%
% INPUTs (all required)
% x     - [1 n] vector vector of, e.g. m/z values
% y     - [1 n] vector of data, e.g. mean spectrum
% col   - [1 n] vector of colours, e.g. PC1 loadings

% Ensure that inputs are [1 n] row vectors...
if size(x,1) > size(x,2)
    x = x';
end
if size(y,1) > size(y,2)
    y = y';
end
if size(col,1) > size(col,2)
    col = col';
end

% Check that matrices are equally sized
if numel(x) ~= numel(y) && numel(x) ~= numel(col)
    error('Matrices are not the same sizes');
end


% Create a new figure
%fig = figure('Position',[100 100 1200 400]); hold on;

% Empty matrix of zeros for z...
z = zeros(size(x));

% Plot the figure using the surface command
surface([x;x],[y;y],[z;z],[col;col],...
    'facecol','no',...
    'edgecol','interp',...
    'linew',2);

ax = gca;
set(ax,'FontSize',14);

% Set the axes limits
xlim([min(x) max(x)]);

box on;

cb = colorbar;
set(cb,'FontSize',14);
if nargin == 4
    ylabel(cb,cbLab);
end
%caxis([min(col) max(col)]);

% return
%     
%     
% %set(gca,'XDir','reverse');
% xlabel('m/z','FontSize',14);
% ylabel(yLab, 'FontSize',14);
% colormap(jet);
% cb = colorbar;
% ylabel(cb,cLab,'FontSize',14);
% 
% 
% % Save the image if a file name has been provided...
% if ~isempty(figNam)
%     graphFormat(figNam,'jpeg');
% end

end

