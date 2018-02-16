function [ ax,centr ] = scatEll(fig,sc,gr,comp)
%scatEll - make a nice scatter plot of x,y and groups, with appropriate
%colours and legend etc... If fig is [], then it draws a new one, otherwise
%it will put it in the provided figure...
allSymb = {'^','d','s','h','v','<','>','p'};

% Error check on the 'fig' variable

% Determine the groupings
[unq,~,ind] = unique(gr)
numG = numel(unq);

if isa(gr,'double')
    axLab = cell(numG,1);
    flag = true;
else
    flag = false;
end

% Somewhere to save group centroids
centr = zeros(numG,2);

% Define colours - may want better than this in the future
if numG > numel(allSymb)
    symb = repmat({'o'},[1 numG]);
    cols = jet(numG);
else
    symb = allSymb(1:numG);
    cols = repmat([0 0 0],[numG 1]);
end

% Create a new figure, or place in the existing one...
if isempty(fig)
    figure; hold on;
else
    axes(fig);
    cla reset;
    hold on;
end

% Now we can draw...
h = zeros(numG,1);
for n = 1:numG
    
    % Group indices
    fx = ind == n;
    
    % Determine the ellipse (shall we say 95%?)
    [ell,centr(n,:)] = error_ellipse(sc(fx,:),95);
    
    % Plot the ellipse
    plot(ell(:,1),ell(:,2),'-',...
        'Color',cols(n,:),...
        'LineWidth',2);
    
    % Add a massive centroid
    h(n,1) = scatter(centr(n,1),centr(n,2),150,cols(n,:),symb{n},'filled',...
        'MarkerEdgeColor','black');
    
    % Now how about the mini scatter points?
    scatter(sc(fx,1),sc(fx,2),50,cols(n,:),symb{n});%,'filled',...
        %'MarkerEdgeColor','black');
    
    if flag
        axLab{n,1} = ['Group ' int2str(unq(n))];
    end
end

ax = gca;

% Add a legend
if flag
    legend(h,axLab,'FontSize',16);
else
    legend(h,unq,'FontSize',16);
end
    
% Format the axes properly...
if nargin > 3
    xlabel(['LV' int2str(comp(1))],'FontSize',16,'FontWeight','bold');
    ylabel(['LV' int2str(comp(2))],'FontSize',16,'FontWeight','bold');
end
set(gca,'FontSize',14);

box on;

    



end

