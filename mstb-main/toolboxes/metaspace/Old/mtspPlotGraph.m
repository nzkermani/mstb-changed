function mtspPlotGraph(gr,group,xyData)
%mtspPlotGraph - plot connectivity graph...

figure('Position',[423 629 920 574]);

% Plot with a circular layout - subject to change with more nodes
h = plot(gr,'Layout','force','NodeLabel',[]);

h.MarkerSize = 10;

% If ydata is specified, then we just replace parts of h
if nargin == 3
    h.XData = xyData(:,1)';
    h.YData = xyData(:,2)';
end

% Colour the nodes
[h] = colourNodes(h,gr.Nodes.(group));

% Colour the edges
[h] = colourEdges(gr,h);

% Format the axes
set(gca,'Units','normalized',...
    'Position',[0 0.05 0.825 0.9]);
axis off;





end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = colourNodes(h,group)

% Get unique members of this group
[unq,~,ind] = unique(group);

% Update the node colour here
h.NodeCData = ind;

% Add a colorbar for visibility purposes
cb = colorbar;
set(cb,...
    'YTick',1:numel(unq),...
    'YTickLabel',unq,...
    'FontSize',16,...
    'FontWeight','bold',...
    'YDir','reverse');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = colourEdges(gr,h)

% Determine line widths
lw = gr.Edges.Weight;
lw = 5 * (lw / max(lw)) .^ 4;
h.LineWidth = lw;

% Set colours of edges to be white if low and black if high
cw = gr.Edges.Weight;
cw = 90 * (cw / max(cw)) .^ 4;
cw = round(cw) + 1;
cmap = flipud(gray(101));
newc = cmap(cw,:);

% Update here
h.EdgeColor = newc;
h.EdgeAlpha = 0.75;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function legacyCode

% Marker colours to match m/z...
varNames = gr.Nodes.Properties.VariableNames;
if strcmp(varNames{2},'mz')
    flag = 'features';
    
    
    %     mz = cell2mat(gr.Nodes.mz);
    %     h.NodeCData = mz;
    %     cb = colorbar;
    %     set(cb,...
    %         'FontSize',16,...
    %         'FontWeight','bold');
    %     ylabel(cb,'m/z','FontSize',16,'FontWeight','bold');
    
    
    [unq,~,ind] = unique(gr.Nodes.Class);
    h.NodeCData = ind;
    cb = colorbar;
    set(cb,...
        'YTick',1:numel(unq),...
        'YTickLabel',unq,...
        'FontSize',16,...
        'FontWeight','bold')
    
elseif strcmp(varNames{2},'Tissue')
    flag = 'tissue';
    
    [unq,~,ind] = unique(gr.Nodes.Tissue);
    if numel(unq) == 1 && any(strcmp(varNames,'Charge'))
        [unq,~,ind] = unique(gr.Nodes.Charge);
    end
    
    h.NodeCData = ind;
    
    cb = colorbar;
    set(cb,...
        'YTick',1:numel(unq),...
        'YTickLabel',unq,...
        'FontSize',16,...
        'FontWeight','bold')
    
end

set(gca,'Units','normalized',...
    'Position',[0 0.05 0.825 0.9]);
axis off;

% Determine line widths
lw = gr.Edges.Weight;
lw = 5 * (lw / max(lw)) .^ 4;
h.LineWidth = lw;

% Set colours of edges to be white if low and black if high
cw = gr.Edges.Weight;
cw = 90 * (cw / max(cw)) .^ 4;
cw = round(cw) + 1;
cmap = flipud(gray(101));
newc = cmap(cw,:);
h.EdgeColor = newc;

h.EdgeAlpha = 0.75;

%return

% Add node labels?
switch flag
    
    case 'features'
        numN = size(gr.Nodes,1);
        if numN < 5000
            for n = 1:numN
                
                lab = gr.Nodes.Chemical{n};
                if ~strcmp(lab,'Unknown')
                    text(h.XData(n),...
                        h.YData(n),...
                        ['  ' gr.Nodes.Chemical{n}],...
                        'FontSize',20,...
                        'Rotation',5,...
                        'Color','red');
                end
                
            end
        end
        
    case 'tissue'
        
end



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%