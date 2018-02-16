function dpnCorrelationFigure(crr)
%dpnCorrelationFigure - plot a figure with all of the correlation images on
%it...

% Names and number of interpolation methods?
methods = fieldnames(crr)
numM = numel(methods)

% Determine min/max for caxis colour limits
lohi = NaN(numM,2);
for n = 1:numM
    
    tmp = crr.(methods{n});
    lohi(n,:) = [min(tmp(:)) max(tmp(:))];
    
end
ll = [min(lohi(:,1)) max(lohi(:,2))];

% How about somewhere to save all of the reshapen data
numE = prod(size(crr.(methods{1})));
resh = NaN(numE,numM);

% Histograms
xvals = 0:0.01:1;
hsts = NaN(numM,numel(xvals));

% Plot a figure
figure;
ax = zeros(numM,1);
for n = 1:numM
    
    ax(n,1) = subplot(1,numM,n);
    
    % Get the data...
    tmp = crr.(methods{n});
    
    % Reshape for boxplots
    resh(:,n) = tmp(:);
    
    % Determine histogram
    hsts(n,:) = hist(tmp(:),xvals);
    
    % Set values less than 0.8 to be NaN
    mask = tmp < 0.8;
    tmp(mask) = NaN;    
    
    imagesc(tmp);
    title(methods{n});
    caxis([0.8 1])
    
end
    
    



end

