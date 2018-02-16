function makeConfMatPredict(cm,txtAct,txtPred)
% Make a uniform looking confusion matrix

% Calculate percentages
prc = 100 * bsxfun(@rdivide,cm,sum(cm,2));

% Size of the matrices
[sz1,sz2] = size(prc);

% General colour map
cmap = flipud(gray(100));

if sz1 >= 5
    fs = 18;
else
    fs = 24;
end

figure; hold on;

imagesc(prc);
colormap(cmap);
caxis([0 100]);

set(gca,'YDir','reverse',...
    'XAxisLocation','top',...
    'XTick',1:sz2,...
    'XTickLabels',txtPred,...
    'YTick',1:sz1,...
    'YTickLabels',txtAct,...
    'FontSize',18);
    


for n = 1:sz2
    
    for r = 1:sz1
        
        % Skip empty pixels
        if cm(r,n) == 0
            continue;
        end
        
        
        if prc(r,n) > 60
            tc = 'white';
        else
            tc = 'black';
        end
        
        tl = [sprintf('%0.1f',prc(r,n)) '%' char(10) sprintf('%d', cm(r,n))];
        
        text(n,r,tl,...
            'Color',tc,...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...            %'FontUnits','pixels',...
            'FontSize',fs,...
            'FontWeight','bold');
        
    end
    
end

% Add unclassified pixels here
if size(prc,2) == sz1 + 1
    for r = 1:sz1
        % Skip empty pixels
        if qty(r,end) == 0
            continue;
        end
                
        if prc(r,end) > 60
            tc = 'white';
        else
            tc = 'black';
        end
        
        tl = [sprintf('%0.1f',prc(r,end)) '%' char(10) sprintf('%d', qty(r,end))];
        
        text(sz+1,r,tl,...
            'Color',tc,...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...            %'FontUnits','pixels',...
            'FontSize',24,...
            'FontWeight','bold');
        
    end
end        
    
% Axis formatting stuffhere
ylabel('Actual Class   ','FontSize',20,'FontWeight','bold');
title('Predicted Class   ','FontSize',20,'FontWeight','bold');
box on;
axis tight square
set(gca,'LineWidth',5,...
    'TickLength',[0 0]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
