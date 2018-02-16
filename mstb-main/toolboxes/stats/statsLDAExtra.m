function statsLDAExtra(~,~,fig,window)
% statsPCAExtra - plot some things in the extra axes

% Get gui data
sts = guidata(fig.fig);
val = window.extraPlot.Value;
str = window.extraPlot.String;
type = str{val};

% Determine if we actually have a confusion matrix to plot
if isempty(sts.res.mmc.cm)
    type = 'none';
end

% Define the parent axes
parent = fig.ax.conf(1);

% Get the desired font and marker sizes
opts.fS = str2double(window.fontSize.String);
opts.mS = str2double(window.markerSize.String);

% Delete existing objects
f0 = get(parent,'Children');
delete(f0);

% Hold the axes
axes(parent);
hold on;

% Switch through the options...
switch type
    case 'Confusion'
        statsPlotConfMat(sts.res.mmc,parent,opts);
        
        %perfcurve(sts.res.mmc.grp,sts.res.mmc.scores,'Adenoma')
        
    otherwise
        return        
end

set(gca,'XTickLabel',[],...
    'YTickLabel',[],...
    'XTick',0,...
    'YTick',0,...
    'LineWidth',5,...
    'TickLength',[0 0]);

grid off;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotConfMat(res,parent,opts)
% Plot the confusion matrix in the correct axes - this will require a bit
% of good coding to ensure that we can use the axes for something other
% than an image afterwards...

% Define the colour map
cmap = flipud(gray(100));


% Get the conf mat
cm = res.cm.cm;
sz = size(cm,1);
prc = bsxfun(@rdivide,cm,sum(cm,2)) * 100;

% Create a conf mat based on the colour scheme...
img = ones(sz,sz,3);
for n = 1:sz
    for r = 1:sz
        if round(prc(n,r)) > 0
            img(n,r,:) = cmap(round(prc(n,r)),:);
        end
    end
end

% Here draw the confusion matrix
cmax = imagesc(parent,img);

% Set the axes limits
xlim([0.48 sz+0.52]);
ylim([0.48 sz+0.52]);

% Colour map
colormap(parula);
caxis([0 100]);

% Reverse y direction
set(gca,'YDir','reverse',...
    'XTick',[],...
    'YTick',[]);

% Get the colours for the observations...
names = res.cm.names;
numG = numel(names);
cols = parula(numG);

hold on;

% So now let's loop through and place the labels...
for n = 1:sz    
    for r = 1:sz
        
        % On the diagonal we need to label the box with colour/symbol and
        % text for group name
        if r == n && opts.mS > 0
            
            scatter(r,n,opts.mS * 1.5,'k','o','filled');
            
            scatter(r,n,opts.mS,cols(n,:),'o','filled',...
                'MarkerEdgeColor','w');
        end
                
        % Skip empty pixels
        if cm(r,n) == 0
            continue;
        end
                
        % Only do if we have a large font size
        if opts.fS > 0

            % Change the font colour
            if prc(r,n) > 60
                tc = 'white';
            else
                tc = 'black';
            end

            
            % Prepare the text string        
            tl = [sprintf('%0.1f',prc(r,n)) '%' char(10) ...
                char(10) sprintf('%d', cm(r,n))];

            % Place the text
            text(n,r,tl,...
                'Color',tc,...
                'HorizontalAlignment','center',...
                'VerticalAlignment','middle',... 
                'FontSize',opts.fS,...
                'FontWeight','bold');
        end
        
    end
    
end
    
    



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotLoads(mz,ll)


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
