function [ output_args ] = makeClassificationMatrix(prc,txt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

figure; hold on;

img = imagesc(prc)

cmap = flipud(gray(100));
colormap(cmap);
caxis([min(prc(:)) max(prc(:))]);

% Reverse
cax = gca;
set(cax,'YDir','reverse');

% Add text
for n = 1:size(prc,1)
    
    for r = 1:size(prc,2)
        
        % What should be displayed?
        if n == r
            lab = txt{n};
            col = 'white';
            
            patch([0.1 0.9 0.9 0.1]+n-0.5,[0.1 0.1 0.9 0.9]+n-0.5,...
                'blue','FaceAlpha',0.5,...
                'EdgeColor','blue');
            
        else
            
            if prc(n,r) > 50
                col = 'white';
            else
                col = 'black';
            end
            lab = int2str(prc(n,r));
                        
        end
        
        text(r,n,lab,...
            'FontSize',18,...
            'HorizontalAlignment','center',...
            'Color',col);
        
    end
    
end



% Axis formatting stuffhere
%ylabel('Actual Class   ','FontSize',20,'FontWeight','bold');
%title('Predicted Class   ','FontSize',20,'FontWeight','bold');
box on;
axis tight square
set(cax,'LineWidth',5,...
    'TickLength',[0 0],...
    'XTick',[],...
    'YTick',[]);


end

