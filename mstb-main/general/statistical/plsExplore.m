function plsExplore(mz,data,meta,pls,comps)
%plsExplore - run through the PLS results showing scores, loadings and box
%plots of selected variables... Exciting!

if nargin == 4
    comps = [1 2];
end

if numel(comps) == 1
    c1 = comps;
    if c1 > 1
        comps = [1 c1];
    else
        comps = [c1 c1+1];
    end
end



% First let's plot a figure of scores and loadings...
figure('Units','normalized','Position',[0.2 0.6 0.6 0.4]);

% Plot the scores in subplot 1
plotScores2(mz,data,meta,pls,comps)

% PLot the loadings in subplot 2
plotLoadings(mz,data,meta,pls,comps);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotScores(mz,data,meta,pls,comps)

ax1 = subplot(1,2,1); hold on;
grp = sum(pls.comps == 1,2) == 1;
[unq,~,ind] = unique(meta.Genus(grp));
numG = numel(unq);
h = zeros(numG,1);
cols = jet(numG);
for n = 1:numG
    
    fx = ind == n;
    h(n,1) = scatter(pls.xs(fx,comps(1)),pls.xs(fx,comps(2)),80,cols(n,:),...
        'o','filled','MarkerEdgeColor','black');
    
end
legend(h,unq);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotScores2(mz,data,meta,pls,comps)
ax1 = subplot(1,2,1); hold on;

% Define symbols
allSymb = {'^','s','d','v','h'};

% Determine the unmixed samples
grp = sum(pls.comps == 1,2) == 1;
[unq,~,ind] = unique(meta.Genus(grp));
numG = numel(unq);
pdata = pls.scPred(grp,:);

% Now just draw the pure samples
leg1 = zeros(numG,1);
for n = 1:numG
    
    % Find those indices
    fx = ind == n;
    
    leg1(n,:) = scatter(pdata(fx,comps(1)),pdata(fx,comps(2)),...
        200,'black',allSymb{n},'filled','MarkerEdgeColor','white');
    
end

% Now let's shove everything else on...
[unq2,~,ind2] = unique(meta.Genus(~grp));
numU = numel(unq2);
pdata = pls.scPred(~grp,:);

% Define the colours
cols = jet(numU);

leg2 = zeros(numU,1);
for n = 1:numU
    
    fx = ind2 == n;
    
    leg2(n,:) = scatter(pdata(fx,comps(1)),pdata(fx,comps(2)),...
        80,cols(n,:),'o','filled','MarkerEdgeColor','k');
end

legend([leg1; leg2],[unq; unq2]);%,'Location','EastOutside')
    
box on;
%grid('on');

%set(gca,'FontSize',16);

return

h = zeros(numG,1);
cols = jet(numG);
for n = 1:numG
    
    fx = ind == n;
    h(n,1) = scatter(pls.xs(fx,comps(1)),pls.xs(fx,comps(2)),80,cols(n,:),...
        'o','filled','MarkerEdgeColor','black');
    
end
legend(h,unq);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotLoadings(mz,data,meta,pls,comps)
% Plot the loadings in a second axes

ax2 = subplot(1,2,2); hold on;

scatter(pls.xl(:,comps(1)),pls.xl(:,comps(2)),80,mz,'o','filled',...
    'MarkerEdgeColor','black');

title('Click next to (not on) a point to display a boxplot')

% This is what happens when we click near a loading value
set(gca,'ButtonDownFcn',{@doBP,mz,data,meta,pls,comps});

cb = colorbar;
ylabel(cb,'m/z ','FontSize',16,'FontWeight','bold','FontAngle','italic');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doBP(src,event,mz,data,meta,pls,comps)

% Get the current point...
click = get(src,'CurrentPoint');

x = click(1,1);
y = click(1,2);

% Now need to find the loading that is neareset...
dst = bsxfun(@minus,pls.xl(:,comps),[x y]) .^ 2;
dst = bsxfun(@rdivide,dst,max(dst,[],1));

% Find the nearest...
euc = sqrt(dst(:,1) + dst(:,2));
[~,idx] = min(euc);

% Find existing boxplot figures and use them...
f0 = findobj('Tag','plsExpBP');
if numel(f0) == 0
    figure('Tag','plsExpBP','Units','normalized','Position',[0.2 0.1 0.6 0.4]);
else
    figure(f0);
    clf reset;
    set(gcf,'Tag','plsExpBP');

end

% Now how about some kind of ANOVA? Just between the 'pure' samples
sel = sum(pls.comps == 1,2) == 1;
[p,~,stats] = anova1(data(sel,idx),meta.Genus(sel),'off');
sqp = false(size(stats.gnames,1));
pthresh = 0.001;
if p < pthresh
    [hsd] = multcompare(stats,'display','off');
    
    for n = 1:size(hsd,1)
        sqp(hsd(n,1),hsd(n,2)) = hsd(n,end) < pthresh;
        sqp(hsd(n,2),hsd(n,1)) = hsd(n,end) < pthresh;
    end    
    
end

% Determine the statistics
[medVal,gnames] = grpstats(data(:,idx),meta.Genus,{'median','gname'});
[~,srt] = sort(medVal);
grpOrd = cell(numel(srt),1);
for n = 1:numel(srt)    
    grpOrd(n,1) = gnames(srt(n),1);
    
end
    
% Draw the boxplot
boxplot(data(:,idx),meta.Genus,...
    'Orientation','horizontal',...
    'GroupOrder',grpOrd);

% Labels
xlabel('Intensity','FontSize',14,'FontWeight','bold');
titText = ['m/z = ' sprintf('%0.3f',mz(idx)) ', ANOVA p = ' sprintf('%0.3e',p)];
title(titText,'FontSize',18,'FontWeight','bold');


% Add the significance
f2 = axes('Position',[0.75 0.15 0.2 0.2]);
imagesc(sqp);
axis square
cmap = colormap(gray);
colormap(flipud(cmap));
set(gca,'YTickLabel',stats.gnames,'YTick',1:numel(stats.gnames));
set(gca,'XTick',[],'FontSize',14);
title(['p < ' sprintf('%0.3f',pthresh)],'FontWeight','bold');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%