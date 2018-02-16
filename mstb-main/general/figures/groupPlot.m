function groupPlot(x,y,g,type)
%groupPlot - colour observations in y according to the groupings in g

% Can modify to plot just averages as well
if nargin == 3
    type = 'normal';
end

[x,y] = insertZeros(x,y,0.01);

% Determine unique groupings
[unqG,~,unqI] = unique(g);

% How many unique groups are there?
numG = numel(unqG);

% Define standard colour maps
cols = jet(numG);
if numG == 3
    cols = [1 0 0; 0 1 0; 0 0 1];
end

% Pre-allocate handle vector & labels
hand = zeros(numG,1);
labs = cell(numG,1);

% Create an empty figure
figure; hold on;

% Loop through each group
for n = 1:numG
    
    % Logical indices of group n
    fx = unqI == n;
    
    % Now plot the spectra...
    switch type
        
        case {'mean','average'}
            ttt = nanmean(y(fx,:),1);
            tmp = plot(x,ttt,'Color',cols(n,:),'LineWidth',1);
            
        case 'median'
            ttt = nanmedian(y(fx,:),1);
            tmp = plot(x,ttt,'Color',cols(n,:),'LineWidth',1);
        
        otherwise
            tmp = plot(x,y(fx,:),'Color',cols(n,:),'LineWidth',1);
            
    end
        
    % Save the handle for the plot
    hand(n,1) = tmp(1);
    
    % Group label - could just use unqG...
    labs{n,1} = unqG{n};
    
end

% Add the legend
legend(hand,labs);

% Put the box around the axes
box on;

% Add x and y axes labels, font size of 16

% Enlarge axes label font size to 14

end

