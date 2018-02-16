function jsmBoxPlot(y,grp,varargin)
%jsmBoxPlot - see below for inputs and optional parameters
%
%
% INPUTs
% y     [m x n] matrix of data
% grp   {m x 1} cell of histological labels, i.e. different boxes
%
% OPTIONAL
% cols  [g x 3] rgb values of colours for correct number of unique(grp)
% order [g x 1] vector to change order of appearance
% mult  [g x 1] vector to which to multiply the variables


% Determine optional input arguments...
[opts] = getVarArgin(grp,varargin);

% Determine groups...
[unq,~,ind] = unique(grp);
numG = numel(unq);

% Determine the median / percentiles of the data.
vals = NaN(numG,3);
for n = 1:numG
    fx = ind == n;
    vals(n,:) = prctile(y(fx),[25 50 75]);
end


% Override sort order if provided...
if ~isempty(opts.order)
    srtIdx = opts.order;
else
    % Sort...
    [~,srtIdx] = sort(vals(:,2));
    opts.labels = opts.labels(srtIdx);
end

% Some place for the legend
legHand = zeros(numG,1);

% Draw the figure
switch opts.orient
    case 'horizontal'
        figure('Position',[0 0 700 1400]); hold on;
    case 'vertical'
        figure('Position',[0 0 700 350]); hold on;
    otherwise
        axes(opts.orient);
        hold on;
        opts.orient = 'vertical';
end

% Loop through the groups
for n = 1:numG
    
    % Which group to focus on...
    i = srtIdx(n);
    
    % Group index
    fx = ind == i;
    
    % Patch coordinates to contain the data...
    px = [vals(i,1) vals(i,3) vals(i,3) vals(i,1)];
    py = [n-0.4 n-0.4 n+0.4 n+0.4];
        
    % Random indices for the scatter points...
    inds = ((rand(sum(fx),1) - 0.5) * 0.50) + n;
    %inds = n + linspace(-0.4,0.4,sum(fx));
    %[s1,s2] = sort(y(fx),'descend');
    %inds2 = inds(s2);
    
    % Draw the bits and pieces
    switch opts.orient
        case 'horizontal'
            patch(px,py,[0.8 0.8 0.8],'EdgeColor','none');
            
            line([vals(i,2) vals(i,2)],[py(1) py(3)],...
                'Color','k','LineWidth',3);    
            
            legHand(i,1) = scatter(y(fx),inds,60,opts.cols(i,:),...
                'o','filled',...
                'MarkerEdgeColor','black');

            
        case 'vertical'
            patch(py,px,[0.8 0.8 0.8],'EdgeColor','none');
            
            line([py(1) py(3)],[vals(i,2) vals(i,2)],...
                'Color','k','LineWidth',3);
            
            legHand(i,1) = scatter(inds,y(fx),60,opts.cols(i,:),...
                'o','filled',...
                'MarkerEdgeColor','black');

    end
    
    
    % How about whiskers?
    %w1 = vals(n,1) - 1.5 * (vals(n,3)-vals(n,1));
    %w3 = vals(n,3) + 1.5 * (vals(n,3)-vals(n,1));
    %line([w1 vals(n,1)],[n n],'Color','k','LineWidth',3);
    %line([w3 vals(n,3)],[n n],'Color','k','LineWidth',3);
    
        
    
    
end

% Axes properties
set(gca,'FontSize',14,'FontName','Calibri');
box on;

% What are the default axes limits
lowVal = min(y(:));
if lowVal == 0
    lowVal = -0.025 * max(y(:));
end
highVal = max(y(:)) * 1.025;

% Labels and stuff
switch opts.orient
    case 'horizontal'
        
        if ~opts.legend
            set(gca,'YTick',1:numG,...
                'YTickLabel',opts.labels,...
                'FontSize',18,...
                'FontWeight','bold');
        else
            set(gca,'YTick',[]);
        end
        
        set(gca,'YDir','reverse');        
        ylim([0.5 numG+0.5]);
        xlim([lowVal highVal]);

        xlabel(opts.yLabel,...
            'FontWeight','bold',...
            'FontSize',18,...
            'FontName','Calibri');


    case 'vertical'
        
        if ~opts.legend
            set(gca,'XTick',1:numG,...
                'XTickLabel',opts.labels,...
                'FontSize',18,...
                'FontWeight','bold');
        else
            set(gca,'XTick',[]);
        end
        %set(gca,'YDir','reverse');        
        xlim([0.5 numG+0.5]);
        ylim([lowVal highVal]);

        ylabel(opts.yLabel,...
            'FontWeight','bold',...
            'FontSize',18,...
            'FontName','Calibri');

end

% Add a title?
if ~isempty(opts.title)
    title(opts.title,...
        'FontSize',20,...
        'FontWeight','bold');
end


% Are we actually going to be adding a legend?
if opts.legend
    
    % Reorder (even if placebo sort)
    legHand = legHand(srtIdx);
    
    % Add the legend, outside is just plain easiest
    legend(legHand,opts.labels,...
        'Location','SouthOutside',...
        'Orientation','horizontal',...
        'FontSize',18,...
        'FontWeight','bold');
    
    legend boxoff;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = getVarArgin(grp,argsin)
% Defaults for the processing parameters...

% Define the defaults here...
opts.cols       = parula(numel(unique(grp)));
opts.order      = [];
opts.orient     = 'horizontal';
opts.legend     = false;
opts.numDP      = 3;
opts.labels     = [];
opts.yLabel     = 'Arbitrary Intensity';
opts.title      = [];

for i = 1:2:numel(argsin)
    if strcmpi('cols',argsin{i}) ...
            || strcmpi('colours',argsin{i}) ...
            || strcmpi('colour',argsin{i}) ...
            || strcmpi('col',argsin{i})
        opts.cols = argsin{i+1};
        
    elseif strcmpi('order',argsin{i})
        opts.order = argsin{i+1};
        
    elseif strcmpi('orientation',argsin{i}) ...
            || strcmpi('orient',argsin{i})        
        opts.orient = argsin{i+1};
        
    elseif strcmpi('muliplier',argsin{i}) ...
            || strcmpi('mult',argsin{i})
        opts.mult = argsin{i+1};
    
    elseif strcmpi('legend',argsin{i})
        opts.legend = argsin{i+1};

    elseif strcmpi('labels',argsin{i}) ...
            || strcmpi('label',argsin{i})
        opts.labels = argsin{i+1};
        
    elseif strcmpi('ylabel',argsin{i})
        opts.yLabel = argsin{i+1};

    elseif strcmpi('title',argsin{i})
        opts.title = argsin{i+1};

        
    end
end

% Orientation validation
switch lower(opts.orient)
    case {'vertical','vert','v'}
        % Thens this is fine
        opts.orient = 'vertical';
        
    case {'horizontal','hori','h'}
        opts.orient = 'horizontal';
    
end

% Create labels if non existent
if isempty(opts.labels)
    opts.labels = cell(size(opts.cols,1),1);
    for n = 1:size(opts.cols,1)
        opts.labels(n,1) = {['Group ' int2str(n)]};
    end
end
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
