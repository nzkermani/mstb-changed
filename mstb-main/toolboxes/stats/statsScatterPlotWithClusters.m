function [patchHand] = statsScatterPlotWithClusters(parent,scores,histID,doEllipse,doCentroids,meta,clust)

% SCatter plot the MVA results

maxGroups = 200;

% Delete existing objects
f0 = get(parent,'Children');
delete(f0);
legend(parent,'off');

% Hold the axes
axes(parent);
hold on;

% Determine unique groups in the histIDs
[unq,~,ind] = unique(histID);
numG = numel(unq);
cols = parula(numG);
symb = 'o';
grpInfo = cell(numG,2);
axLeg = zeros(numG,1);

if isnumeric(unq)
    unq = num2str(unq);
end
    
if numG <= maxGroups
    for n = 1:numG

        % Indices of points to plot
        fx = ind == n;

        % Colour finding/matching
        col = cols(n,:);

        grpInfo{n,1} = col;
        grpInfo{n,2} = symb;
        
        % Determine centroid (whether or not we need it)
        cntrd = mean(scores(fx,:),1);

        % Scatter away
        if ~doCentroids
            axLeg(n,1) = scatter(scores(fx,1),scores(fx,2),...
                50,...
                col,...
                symb,...
                'filled',...
                'MarkerEdgeColor','k',...
                'HitTest','on',...
                'ButtonDownFcn',{@quickScatterClick,parent,scores(:,1),scores(:,2),meta});

        else
            axLeg(n,1) = scatter(cntrd(1),cntrd(2),...
                120,...
                col,...
                symb,...
                'filled',...
                'MarkerEdgeColor','k');
        end
            

        % Determine the ellipse
        if numG < maxGroups && sum(fx) >= 5 && doEllipse
            [ell,~,~] = error_ellipse(scores(fx,:),95);
            plot(ell(:,1),ell(:,2),'Color',col,'LineWidth',2);
        end
                
    end
    
else
    scatter(scores(:,1),scores(:,2),120,ind,symb,'filled',...
            'MarkerEdgeColor','k',...
                'HitTest','on',...
                'ButtonDownFcn',{@quickScatterClick,parent,scores(:,1),scores(:,2),meta});
end
    
    
    
if numG <= maxGroups
    lgnd = legend(axLeg,unq);
    set(lgnd,'FontSize',14);
else
    legend('hide');
end

% Draw patches
[patchHand] = patchyMcPatchFace(scores,clust);

box on;

set(gca,'XTickLabel',[],'YTickLabel',[],...
    'XTick',0,...
    'YTick',0,...
    'LineWidth',5,...
    'TickLength',[0 0]);

grid off;

% Draw lines along the origin
xlim('auto');
ylim('auto');
xl = xlim(gca);
yl = ylim(gca);
line([xl(1) xl(2)],[0 0],'LineStyle',':','Color','k');
line([0 0],[yl(1) yl(2)],'LineStyle',':','Color','k');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quickScatterClick(~,event,parent,xx,yy,histID)

% Remove other references that have been clicked
f0 = findobj('Tag','quickClick');
delete(f0);

% Only show with left click
if event.Button > 1
    return
end

% Get coordinates
coor = get(parent,'CurrentPoint');
x = coor(1,1);
y = coor(1,2);
xd = (xx - x) .^ 2; 
xd = xd / nanmax(xd);
yd = (yy - y) .^ 2; 
yd = yd / nanmax(yd);
[~,dd] = min(xd + yd);

% Format the metadata into a single string
fn = fieldnames(histID);

% Can we find the best couple of examples? Ideally we have fileID and
% histID, but in other instances we need to look for alternatives...
fileI = strcmp(fn,'fileID') | strcmp(fn,'patientID');
histI = strcmp(fn,'histID') | strcmp(fn,'tissueID');

% Combine any matches of these
idx = find(fileI | histI);
if numel(idx) == 0
    idx = [1 2];
end

% Prepare the string, beware as some might be numeric
mdstr = histID.(fn{idx(1)}){dd};
for n = 2:numel(idx)
    mdstr = [mdstr ', ' histID.(fn{idx(n)}){dd}];
end

% Now add in the text...
span = xlim;
span = span(2) - span(1);
xPos = xx(dd) + span * 0.01;
text(xPos,yy(dd),mdstr,'Tag','quickClick',...
    'FontSize',14,...
    'FontWeight','bold',...
    'BackgroundColor','green',...
    'Interpreter','none');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [patchHand] = patchyMcPatchFace(scores,clust)

% Unique clusters
[unq,~,ind] = unique(clust);
numC = numel(unq);
patchCols = parula(numC);
patchHand = zeros(numC,1);

% Determine outline for each cluster and then draw it as a patch on the
% figure...
for n = 1:numC
    
    % Indices for cluster n
    fx = ind == n;
    
    % Points in this cluster
    pxy = scores(fx,1:2);
    
    % Cluster boundary
    if size(pxy,1) == 1
        pcx = pxy(1) + [-10 10 10 -10];
        pcy = pxy(2) + [-10 -10 10 10];
    else
        ol = boundary(pxy(:,1),pxy(:,2));
    
        % Patch coordinates
        pcx = pxy(ol,1);
        pcy = pxy(ol,2);
    end
    
    patchHand(n,1) = patch(pcx,pcy,'red',...
        'FaceColor',patchCols(n,:),...
        'FaceAlpha',0.2,...
        'Tag',['Patch' int2str(n)]);
    
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
