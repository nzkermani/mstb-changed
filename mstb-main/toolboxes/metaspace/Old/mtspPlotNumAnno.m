function mtspPlotNumAnno(numAnno,grp)
%mtspPlotNumAnno

[unq,~,ind] = unique(grp);
numG = numel(unq);
unq

figure;

% Histogram limits...
lims = [0 max(numAnno)];
step = 10;
posn = [0:step:(lims(2)+step)];

ax = zeros(numG,1);

% Loop through...
for n = 1:numG
    
    fx = ind == n;
    
    [y,x] = hist(numAnno(fx),posn);
    
    ax(n,1) = subplot(1,numG,n);
    bar(x,y);
    
    xlim([-1 posn(end)]);
    
    set(gca,'FontSize',14,'FontWeight','bold');
    title(unq{n},'FontSize',16,'FontWeight','bold');
    
    xlim(ax(n,1),[-1 100]);
end


end

