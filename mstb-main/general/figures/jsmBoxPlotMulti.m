function jsmBoxPlotMulti(y,grp,grp2,varargin)%cols,order,orient)
%jsmBoxPlotMulti - make Multiple box plots. grp is the histological group,
%for example, and grp2 is the variable name
%
%
% INPUTs
% y     [m x n] matrix of data
% grp   {m x 1} cell of histological labels, i.e. different boxes
% grp2  [1 x n] vector of variable names, e.g. m/z values
%
% OPTIONAL
% cols  [g x 3] rgb values of colours for correct number of unique(grp)
% order [g x 1] vector to change order of appearance
% mult  [g x 1] vector to which to multiply the variables

% Determine optional input arguments...
[opts,grp2] = getVarArgin(grp,grp2,varargin);

% Determine groups...
[unq,~,ind] = unique(grp);
numG = numel(unq);

% if nargin > 2
%     if size(cols,1) < numG
%         warning('colours incorrect')
%         cols = jet(numG);
%     end
% else
%     cols = jet(numG);
% end
% 
% if nargin == 4
%     order = [];
%     orient = 'horizontal';
% end

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

% Here we would multiply the variables by the multiplier if it is provided.
% It is set to 1 by default if not specified, so can just scale accordingly
y = bsxfun(@times,y,opts.mult);

% How many variables are there?
numV = numel(grp2);

for q = 1:numV
    
    % Determine the median / percentiles of the data.
    vals = zeros(numG,3);
    for n = 1:numG
        fx = ind == n & ~isinf(y(:,q));
        vals(n,:) = prctile(y(fx,q),[25 50 75]);
    end

    % Override sort order if provided...
    if ~isempty(opts.order)
        srtIdx = opts.order;
    else
        [~,srtIdx] = sort(vals(:,2));
    end

    % Height of individual box
    st = 0.1;
    fn = 0.9;
    ht = 1/(fn-st);
    sp = linspace(st,fn,numG+1);
    bu = 0.005;
    
    for n = 1:numG

        % Which group to focus on...
        i = srtIdx(n);

        % Group indices
        fx = ind == i;

        % Patch to contain the data...
        px = [vals(i,1) vals(i,3) vals(i,3) vals(i,1)];
        py = [sp(n)+bu sp(n)+bu sp(n+1)-bu sp(n+1)-bu] + 0.5 + q-1;      
        
        % Random indices in the y axis... between py(1) and py(end)
        inds = py(1) + (rand(sum(fx),1) * (py(end)-py(1)));

        % Do the plotting in here
        switch opts.orient
            case 'horizontal';
                
                % Patch
                patch(px,py,opts.cols(i,:),'EdgeColor','none','FaceAlpha',0.5);
                
                % Scatter plot the points...
                legHand(i,1) = scatter(y(fx,q),inds,80,opts.cols(i,:),...
                    'o','filled',...
                    'MarkerEdgeColor','black');
        
                % Median line
                line([vals(i,2) vals(i,2)],[py(1) py(3)],...
                    'Color','k',...
                    'LineWidth',3);

            case 'vertical'
                
                % Patch                
                patch(py,px,opts.cols(i,:),...
                    'EdgeColor','none',...
                    'FaceAlpha',0.5);
                
                % Scatter plot the points...
                legHand(i,1) = scatter(inds,y(fx,q),30,opts.cols(i,:),...
                    'o','filled',...
                    'MarkerEdgeColor','black');
        
                % Median line
                line([py(1) py(3)],[vals(i,2) vals(i,2)],...
                    'Color','k',...
                    'LineWidth',3);

        end


%         % Scatter plot the points...
%         legHand(i,1) = scatter(y(fx,q),inds,30,cols(i,:),'o','filled',...
%             'MarkerEdgeColor','black');
%         
%         % Median line
%         line([vals(i,2) vals(i,2)],[py(1) py(3)],...
%             'Color','k',...
%             'LineWidth',3);

    end

end



% Are we actually going to be adding a legend?
if opts.legend
        
    % % Reorder the legend to match the order...
    % if ~isempty(opts.order)
    %     legHand = legHand(opts.order);
    % end

    % Reorder (even if placebo sort)
    legHand = legHand(srtIdx);
    
    % Add the legend, outside is just plain easiest
    legend(legHand,opts.labels,...
        'Location','NorthOutside',...
        'Orientation','horizontal',...
        'FontSize',18,...
        'FontWeight','bold');
    
    legend boxoff;
end


% Axes properties
set(gca,'FontSize',14,...
    'FontName','Calibri');
box on;

% xlabel
if min(y(:)) < 0
    xtext = 'Log_2 Fold Change';
else
    xtext = 'Arbitrary Intensity';
end


switch opts.orient
    case 'horizontal'
        xlabel(xtext,...
            'FontWeight','bold',...
            'FontSize',18,...
            'FontName','Calibri');

        % ylabel
        set(gca,'YTick',1:numV,...
            'YTickLabel',grp2,...
            'FontSize',18,...
            'FontWeight','bold');

        set(gca,'YDir','reverse');

        ylim([0.5 numV+0.5]);
        xlim([min(y(:)) max(y(:))]);
        
    case 'vertical'
        ylabel(xtext,...
            'FontWeight','bold',...
            'FontSize',18,...
            'FontName','Calibri');

        % ylabel
        set(gca,'XTick',1:numV,...
            'XTickLabel',grp2,...
            'FontSize',18,...
            'FontWeight','bold');

        %set(gca,'YDir','reverse');

        xlim([0.5 numV+0.5]);
        ylim([min(y(:)) max(y(:))]);

        
end


% Finally, if we have used the multiplication factor, then we need to
% display that tastefully on the axes...
for n = 1:numV
    
    % Simple multiplication check / abort function
    if opts.mult(n) == 1
        % Then there is nothing to do
        continue;
    end
        
    % Need to determine the x/y location for the text...
    switch opts.orient
        case 'vertical'
            xloc = n;
            yloc = full(max(y(:)) * 0.925);
        case 'horizontal'
            yloc = n;
            xloc = full(max(y(:)) * 0.925);
    end
    
    ts = [char(9) '\times' char(9) char(9) ' ' int2str(opts.mult(n))];
    ts = [char(215) sprintf('%d',opts.mult(n))]
    
    text(xloc,yloc,ts,...
        'VerticalAlignment','middle',...
        'HorizontalAlignment','center',...
        'FontSize',18,...
        'FontWeight','bold');
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts,new] = getVarArgin(grp,grp2,argsin)
% Defaults for the processing parameters...

% Define the defaults here...
opts.cols       = jet(numel(unique(grp)));
opts.order      = [];
opts.orient     = 'horizontal';
opts.mult       = ones(1,numel(grp2));
opts.legend     = false;
opts.numDP      = 2;
opts.labels     = [];

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
        
    elseif strcmpi('labels',argsin{i}) ...
            || strcmpi('label',argsin{i})
        opts.labels = argsin{i+1};

    
    elseif strcmpi('legend',argsin{i})
        opts.legend = argsin{i+1};

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

% Convert the grp2 to a cell string for fixed decimals
if ~iscell(grp2)
    new = cell(numel(grp2),1);
    for n = 1:numel(grp2)
        new{n,1} = sprintf(['%0.' int2str(opts.numDP) 'f'],grp2(n));        
    end
else
    new = grp2;
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
