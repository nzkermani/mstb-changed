function [fig,ax] = scatterPlot(xy,grp,txt,xl,yl,names)
% scatterPlot - function to plot data

[unq,~,ind] = unique(grp);
numG = numel(unq);

if numG == 2
    cols = [0 1 0; 1 0 0];
else
    cols = parula(numG);
end
%cols = [0 0 1; 0 0 0.5; 0 1 0; 0 0.5 0; 1 0 0; 0.5 0 0];
symb = 'o';%,'d','o','d','o','d'};
axH = zeros(numG,1);

if isa(unq,'double')
    unq2 = cell(numG,1);
    flag = true;
else
    flag = false;
end

[fig] = figure('Position',[100 100 800 600]); hold on;

for n = 1:numG
    
    fx = ind == n;
    
    axH(n,1) = scatter(xy(fx,1),xy(fx,2),80,cols(n,:),...
        symb,'filled',...
        'MarkerEdgeColor','k');
    
    % Draw an ellipse
    try
        [ell] = error_ellipse(xy(fx,1:2),95);
        plot(ell(:,1),ell(:,2),'Color',cols(n,:),'LineWidth',2);
    catch
        disp('not enough obs');
    end
    
    if ~isempty(txt)
        text(xy(fx,1)+0.00751,xy(fx,2)+0,txt(fx,:));
    end    
    
    if flag
        if isempty(names)
            unq2{n,1} = ['Group ' int2str(unq(n,1))];
        else
            unq2{n,1} = names{n,1};
        end
    end
end
    
if flag
    legend(axH,unq2,'Location','NorthEastOutside');
else
    legend(axH,unq,'Location','NorthEastOutside');
end

box on;

set(gca,'FontSize',16);

xlabel(xl,'FontSize',18,'FontWeight','bold');
ylabel(yl,'FontSize',18,'FontWeight','bold');

ax = gca;

end

