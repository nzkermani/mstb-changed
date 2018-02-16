function [fig,ax] = pcaPlot(xy,grp,varargin)
% scatterPlot - function to plot data

% Various arguments
[opts] = readArgsData(varargin);

% Determine the number of groups
[unq,~,ind] = unique(grp);
numG = numel(unq);

% Draw the figure
[fig] = figure('Position',[100 100 800 600]); hold on;

% Decide if this is continuous data, in which case we skip the groups and
% plot everything together.
counts = hist(ind,max(ind));
if (numG == 1 || max(counts) < 3) && isnumeric(grp)
    
    scatter(xy(:,1),xy(:,2),80,grp,...
        'o','filled',...
        'MarkerEdgeColor','k');
    
    if numG == 1
        [ell,~,~] = error_ellipse(xy(:,1:2),opts.ellipseCI); 
        plot(ell(:,1),ell(:,2),'Color','k','LineWidth',2);        
    end

    cb = colorbar;
    
    ylabel(cb,opts.cbLab,'FontSize',opts.fsLab);
    
    % Box is a classic
    box on;

    % Set the font size
    set(gca,'FontSize',opts.fsAxes);

    % Add the labels
    xlabel(opts.xLab,'FontSize',opts.fsLab,'FontWeight','bold');
    ylabel(opts.yLab,'FontSize',opts.fsLab,'FontWeight','bold');

    ax = gca;

    return
end

% Somewhere to store the centroids...
cent = zeros(numG,2);

% This is for numeric entries... as the legend cannot handle numeric
% entries
if isa(unq,'double')
    unq2 = cell(numG,1);
    flag = true;
else
    flag = false;
end

% Define colours...
try
    cols = parula(numG);
catch
    cols = jet(numG);
end

% Axes handles
axH = zeros(numG,1);

% Scatter each group
for n = 1:numG
    
    % Group indices
    fx = ind == n;
        
    % See if we need to determine the ellipse
    if sum(fx) > 2 && (opts.ellipse || opts.centroid)
        [ell,cent(n,:),~] = error_ellipse(xy(fx,1:2),opts.ellipseCI); 
    end
    
    % Plot the ellipse
    if sum(fx) > 2 && opts.ellipse
        plot(ell(:,1),ell(:,2),'Color',cols(n,:),'LineWidth',2);        
    end

    % Scatter the points?
    if opts.points
        axH(n,1) = scatter(xy(fx,1),xy(fx,2),80,cols(n,:),...
            'o','filled',...
            'MarkerEdgeColor','k');
    end
    
    % Plot the centroid
    if opts.centroid
        axH(n,1) = scatter(cent(n,1),cent(n,2),80,cols(n,:),...
            'o','filled',...
            'MarkerEdgeColor','k');
    end
      

    % Add optional text labels
    if ~isempty(opts.textLabels)
        text(xy(fx,1),xy(fx,2)+80,opts.textLabels(fx,:));
    end    
    
    % This is needed in case the labels are numeric (or missing)
    if flag
        if isempty(opts.groups)
            unq2{n,1} = ['Group ' num2str(unq(n,1))];
        else
            unq2{n,1} = opts.groups{n,1};
        end
    end
end

% What about the legend?
if opts.legend
    if flag
        legend(axH,unq2,'Location','Best');
    else
        legend(axH,unq,'Location','Best');
    end
end

% Connect the centroids?
if opts.connect
    plot(cent(:,1),cent(:,2),'LineWidth',1,'Color','k');
end

if flag && opts.showCB
    cb = colorbar;
    caxis([1 numG]);
    ylabel(cb,opts.cbLab,'FontSize',opts.fsLab,'FontWeight','bold');
end

% Box is a classic
box on;

% Set the font size
set(gca,'FontSize',opts.fsAxes);

% Add the labels
xlabel(opts.xLab,'FontSize',opts.fsLab,'FontWeight','bold');
ylabel(opts.yLab,'FontSize',opts.fsLab,'FontWeight','bold');

ax = gca;

return







%cols = [0 0 1; 0 0 0.5; 0 1 0; 0 0.5 0; 1 0 0; 0.5 0 0];
%cols = [0 0 1; 1 0 0];
symb = 'o';%,'d','o','d','o','d'};



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
        text(xy(fx,1),xy(fx,2)+80,txt(fx,:));
    end    
    
    if flag
        if isempty(opts.groups)
            unq2{n,1} = ['Group ' int2str(unq(n,1))];
        else
            unq2{n,1} = names{n,1};
        end
    end
end
    

box on;

set(gca,'FontSize',16);

xlabel(xl,'FontSize',18,'FontWeight','bold');
ylabel(yl,'FontSize',18,'FontWeight','bold');

ax = gca;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = readArgsData(argsin)
% Read the arguments and then the data if it wasn't passed

% Define the defaults here
opts.ellipse = true;
opts.ellipseCI = 95;
opts.groups = [];
opts.legend = true;
opts.xLab = 'x-axis';
opts.yLab = 'y-axis';
opts.fsLab = 18;
opts.fsAxes = 16;
opts.markers = {'o','s','d'};
opts.points = true;
opts.cbLab = '';
opts.textLabels = [];
opts.connect = false;
opts.showCB = false;
opts.centroid = false;

% Run through each pair
nArgs = length(argsin);
for i = 1:2:nArgs
    if strcmpi('ellipse',argsin{i})
        tmp = argsin{i+1};        
        if islogical(tmp)
            opts.ellipse = tmp;
        end
        
    elseif strcmpi('confidence',argsin{i})
        tmp = argsin{i+1};
        if isnumeric(tmp)
            if tmp == 90 || tmp == 95 || tmp == 99 || tmp == 99.9
                opts.fP = tmp;
            end
        end
        
    elseif strcmpi('points',argsin{i})
        tmp = argsin{i+1};
        if islogical(tmp)
            opts.points = tmp;
        end
        
    elseif strcmpi('centroid',argsin{i})
        tmp = argsin{i+1};
        if islogical(tmp)
            opts.centroid = tmp;
        end
        
    elseif strcmpi('groups',argsin{i})
        tmp = argsin{i+1};
        opts.groups = tmp;
        
    elseif strcmpi('legend',argsin{i})
        tmp = argsin{i+1};
        if islogical(tmp)
            opts.legend = tmp;
        end
        
    elseif strcmpi('connect',argsin{i})
        tmp = argsin{i+1};
        if islogical(tmp)
            opts.connect = tmp;
        end
        
    elseif strcmpi('cb',argsin{i})
        tmp = argsin{i+1};
        if islogical(tmp)
            opts.showCB = tmp;
        end
        
    elseif strcmpi('xlabel',argsin{i})
        tmp = argsin{i+1};
        if ischar(tmp)
            opts.xLab = tmp;
        end
    elseif strcmpi('ylabel',argsin{i})
        tmp = argsin{i+1};
        if ischar(tmp)
            opts.yLab = tmp;
        end
        
    elseif strcmpi('cblabel',argsin{i})
        tmp = argsin{i+1};
        if ischar(tmp)
            opts.cbLab = tmp;
        end
        
    elseif strcmpi('sizelabel',argsin{i})
        tmp = argsin{i+1};
        if isnumeric(tmp)
            opts.fsLab = tmp;
        end
        
    elseif strcmpi('sizeaxes',argsin{i})
        tmp = argsin{i+1};
        if isnumeric(tmp)
            opts.fsAxes = tmp;
        end
        
    elseif strcmpi('confidence',argsin{i})
        tmp = argsin{i+1};
        opts.markers = tmp;       
    
    elseif strcmpi('text',argsin{i})
        tmp = argsin{i+1};
        opts.textLabels = tmp;       


    
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
