function [ hFig,hAx ] = pca_ell( sc, gr, cl )
%pca_ell - draw a pca plot with ellipses instead of the scatter points...
%Requires you to have done the PCA before, passing the scores and the
%groups around which you want to form the ellipses.

% Decide about multiple classes here...
if nargin == 2
    cl = repmat({'x'},size(sc,1),1);
end

% What about unique classes?
[cUnq,~,cInd] = unique(cl);

% What about unique groups?
[gUnq,~,gInd] = unique(gr);

% Combine for unique rows of cInd and gInd
comb = [cInd gInd];
[aUnq,~,aInd] = unique(comb,'rows');
   
numP = size(aUnq,1);
    
if numel(gUnq) == 4
    cols = [1 0 0; 0 1 0; 0 0 1; 0.5 0 0.5];    
else
    cols = jet(numel(gUnq));
end




hFig = figure; hold on;

xys = zeros(numP,2);
for n = 1:numP
    
    % Extract this group from the scores
    fx = find(aInd == n);
    
    % Get the ellipse...
    [ell,xys(n,:)] = error_ellipse(sc(fx,:));
    
    ci = aUnq(n,2);
    
    % Draw it on the graph
    plot(ell(:,1),ell(:,2),'-','Color',cols(ci,:));
    
    
end

% And the center points?
symb = {'o','s','d','p','h'};
hAx = zeros(numel(cUnq),1);

% Define the colour matrix, and size appropriately
% if size(aUnq,1) ~= size(cols,1)
%     fac = size(aUnq,1) / size(cols,1)
%     cols = repmat(cols,fac,1);        
% end

%return

for n = 1:numel(cUnq)
    % Each gets a different symbol (hopfully!)
    si = rem(n,numel(symb));
    if si == 0
        si = numel(symb);
    end
        
    % Which ones to plot?
    fx = aUnq(:,1) == n;
    
    try
        hAx(n,1) = scatter(xys(fx,1),xys(fx,2),...
            90,repmat(cols,numel(cUnq),1),symb{si},'filled',...
            'MarkerEdgeColor','black');
    catch
        try
            hAx(n,1) = scatter(xys(fx,1),xys(fx,2),...
                90,cols(fx,:),symb{si},'filled',...
                'MarkerEdgeColor','black');
        catch
            hAx(n,1) = scatter(xys(fx,1),xys(fx,2),...
                90,'black',symb{si},'filled',...
                'MarkerEdgeColor','black');
        end
    end

end

legend(hAx,cUnq,'FontSize',14,'Location','NorthEast');

end

