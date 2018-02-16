function [ output_args ] = volcanoPlot(g)
%volcanoPlot

% Get guidata if none
if nargin == 0
    g = guidata(gcf);
    if isempty(g)
        error('No guidata in this figure');
    end
end

% Check that ANOVA has been performed
if ~isfield(g,'BR')
    error('Do ANOVA first');
end

% If we can find a tumour like class, we'll make that group one. Otherwise,
% use just 1 and 2...
g.groupIds
isT = strfind(lower(g.groupIds),'tumour');
isC = strfind(lower(g.groupIds),'cancer');
cls = ~cellfun(@isempty,[isT isC]);
if any(cls(:)) == 1
    [~,cls1] = max(sum(cls,2));
    if cls1 == 1
        cls2 = 2;
    else
        % Uses 1 by default for the other class
        cls2 = 1;
    end
else
    cls1 = 1;
    cls2 = 2;
end

% These are q-values from ANOVA
qVal = -log10(g.BR.qvals); %pvalues
qVal(isinf(qVal)) = 0;
qVal(isnan(qVal)) = 0;

% Find the two classes that we are interested in...
f1 = g.groupdata == cls1;
f2 = g.groupdata == cls2;

% Means of each group
m1 = nanmedian(g.Sp(f1,:),1);
m2 = nanmedian(g.Sp(f2,:),1);

% Log2 fold changes
lfc = log2(m1./m2);
lfc(isinf(lfc)) = NaN;

% Do some basic housekeeping:
% - set all 0-fold change variables to qVal = NaN;
% - set all 0-qVal variables to fold change = NaN;
noFC = lfc == 0;
noQV = qVal == 0;
lfc(noQV) = 0;
qVal(noFC) = 0;

% Determine the best variables...
%[~,rank1] = sort(abs(lfc),'descend');
%[~,rank2] = sort(qVal,'descend');


fig.fig = figure('Name','Volcano!'); hold on;
guidata(fig.fig,g);

% Draw the subaxes...
fig.ax(2) = subplot(1,2,2);
fig.ax(1) = subplot(1,2,1);

% Draw the main volcano plot
scatter(lfc,qVal,40,g.ppm,'o','filled',...
    'MarkerEdgeColor','black',...
    'HitTest','on',...
    'ButtonDownFcn',{@figClick,fig,lfc,qVal}); 

title(['Fold changes: ' g.groupIds{cls1} ' / ' g.groupIds{cls2}],...
    'FontSize',16);
xlabel('log_2 fold change','FontSize',16);
ylabel('-log_10 q-value','FontSize',16);
cbar = colorbar;
ylabel(cbar,'m/z','FontSize',16);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figClick(src,event,fig,lfc,qval)
% Get figure click... and then...

% Get coordinates
coor = get(fig.ax(1),'CurrentPoint');
x = coor(1,1);
y = coor(1,2);

%return

% Get the actual data...
g = guidata(fig.fig);

% Calculate the nearest to the lfc
dFC = (lfc - x) .^ 2; dFC = dFC / max(dFC);
dQV = (qval - y) .^ 2; dQV = dQV / max(dQV);
[~,dd] = min(dFC + dQV);

% Reset the secondary axes
axes(fig.ax(2));
%cla reset;

% Now do the boxplot of this variable...
boxplot(fig.ax(2),g.X(:,dd),g.groups,'Labels',g.groupIds);
ylabel(fig.ax(2),'Logged intensity','FontSize',16);

txt = ['m/z = ' sprintf('%0.4f',g.ppm(dd)) ...
    ', log2FC = ' sprintf('%0.2f',lfc(dd)) ...
    ', q = ' sprintf('%0.4e',10^-qval(dd))];
title(fig.ax(2),txt,'FontSize',16);

set(fig.ax(1),'ButtonDownFcn',{@figClick,fig,lfc,qval});%,'HitTest','off');

drawnow;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%