function [ output_args ] = sigFeat(x,y,g,rank,qty)
%sigFeat - generate plots of significant features.  This currently just
%calls the boxplot function on an iterative basis.
%
% INPUTs
% x     - [1 x n], vector of mz,ppm values
% y     - [m x n], matrix of spectra
% g     - {m x 1}, cell array or vector of groupings/values
% rank  - [n x 1], variable, e.g. loadings
% qty   - [1 x 1], number telling us how many plots to generate


% Find the indices of the qty largest values of rank
[~,idx] = sort(rank,'descend');
idx = idx(1:qty)'

% Loop through the number of desired figures
for n = 1:qty
    
    figure; hold on;
    
    % Determine the correlation
    crr = corr(y(:,idx(n)),g);
    
    % Determine the line of best fit
    cf = polyfit(y(:,idx(n)),g,1);
    
    lx = polyval(cf,0:1:100);
        
    scatter(g,y(:,idx(n)),80,g,'o','filled',...
        'MarkerEdgeColor',[0.5 0.5 0.5]);
    
    %plot(0:1:100,lx,'k');
    
    title(['m/z = ' sprintf('%0.1f',x(idx(n))) ', r = ' sprintf('%0.2f',crr)],...
        'FontSize',18,'FontWeight','bold');
    
    xlabel('Tumour Percentage',...
        'FontSize',18,'FontWeight','bold');
    
    ylabel('Arbitrary Intensity',...
        'FontSize',18,'FontWeight','bold');
    
    set(gca,'FontSize',16);
    box on;
    
    graphFormat(['/Users/jmckenzi/Dropbox/Imperial/Work4Others/David/FigPred' int2str(n)],'png');
end

end

