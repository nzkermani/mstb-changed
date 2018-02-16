function statsUnivariateClick(~,~,fig,lfc,qval,sp,histID,xy,cols,pqString)
% Get figure click... and then...

% Get coordinates
coor = get(fig(1),'CurrentPoint');
x = coor(1,1);
y = coor(1,2);

% Delete existing hitPoints...
f0 = findobj('Tag','hitPoint');
delete(f0);

% xy coordinates...
if size(xy,2) == 1    
    xy = [xy qval xy];
end
   
% Calculate the nearest to the lfc
dFC = (xy(:,1) - x) .^ 2; 
dFC = dFC / nanmax(dFC);
dQV = (xy(:,2) - y) .^ 2; 
dQV = dQV / nanmax(dQV);
[~,dd] = min(dFC + dQV);


% Add a marker for the clicked point
hitPoint = scatter(xy(dd,1),xy(dd,2),...
    600,[0.8392 0.3765 0.3020],'p','filled',...
    'MarkerEdgeColor',[0.8392 0.3765 0.3020]);
set(hitPoint,'Tag','hitPoint');

% Reset the secondary axes
axes(fig(2));
f0 = get(fig(2),'Children');
delete(f0); legend(fig(2),'off');
%cla reset;

% Determine the unique values
[unq,~,ind] = unique(histID);


bpmethod = 'jsm';
switch bpmethod
    
    case 'jsm'

        % Use my default boxplot function
        jsmBoxPlot(sp(:,dd),ind,...
            'Orientation',fig(2),...
            'Colours',cols,...
            'Legend',false,...
            'Labels',unq,...
            'Order',1:numel(unq));
       
        % Change some of the axes properties
        set(fig(2),...
            'YTickLabel',[],...
            'FontSize',10,...
            'FontWeight','normal',...
            'FontName','Helvetica',...
            'YDir','normal',...
            'YTickMode','auto',...
            'YTickLabelMode','auto');
        
        ylabel(fig(2),'');
        %title(fig(2),'');

    case 'trad'

        % Now do the boxplot of this variable...
        boxplot(fig(2),...
            sp(:,dd),...
            ind,... % groups
            'Labels',unq,...
            'Jitter',0.8);

end

% Determine the y limits...
yl = ylim(fig(2));
span = (yl(2) - yl(1)) * 1.05;
textY = yl(1) + ((yl(2) - yl(1)) * 1.015);
yl(2) = yl(1) + span;
ylim(fig(2),yl);

% The volcano plot shows different parts, and different inputs are
% provided.
if size(xy,2) == 4
    xy = xy(:,3:4);
else
    xy = xy(:,3);
end

% What label do we want to show?
mzS = sprintf('%0.4f',xy(dd,1));
if size(xy,2) == 2
    rtS = [', RT = ' sprintf('%0.4f',xy(dd,2))];
else
    rtS = '';
end
pqS = sprintf('%0.4e',10^-qval(dd));
fcS = sprintf('%0.2f',lfc(dd));

txt = ['m/z = ' mzS rtS ', log2FC = ' fcS ', ' pqString ' = ' pqS];


textX = 0.55;%numel(unq) - 0.5;
text(textX,textY,txt,'FontSize',14);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%