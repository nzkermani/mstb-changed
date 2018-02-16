function [ output_args ] = opassanal(op,ass,guistruct)
%opassanal - do tests on the data for variables which are significant

% Which groups of the data are we interested in?
groups = {'TUM','MUC'};
numG = numel(groups);
numO = size(op.XPeaks,1);

% Extract only those indices
idx = zeros(numO,numG);
for n = 1:numG    
    idx(:,n) = strcmp(op.histID,groups{n});
end
idx = sum(idx,2) >= 1;

% Additionally, have to exclude some samples for whatever reason
exclude = {'JLA011_normal_S3-.h5','JLA027_tumour_S4-.h5'};
exc = zeros(numO,numel(exclude));
for n = 1:numG
    exc(:,n) = strcmp(op.sampleID,exclude{n});
end
exc = sum(exc,2) >= 1;

% Final list for inclusion
idx = idx & ~exc

% Now do multivariate analysis on those
if isempty(guistruct)
    
    uiDimensionReduction(op.cmz,...
        full(op.XPeaks(idx,:)),...
        full(op.XPeaksLog(idx,:)),...
        op.histID(idx,:),...
        'sampleIDs',op.sampleID(idx,:),...
        'replicates',op.pID(idx,:));
    
    return
    
    % Now we need the multivariate analysis and ANOVA q-correction to be
    % performed and returned to this function...
    
end

% Now we re-run the function providing the output from the stats toolbox.
% We want to find the signficant loadings from multivariate analysis and
% low q-values from FDR-corrected ANOVA

% Group info
gnam = guistruct.groupIds;
gidx = guistruct.groupdata;

% Q values
qq = guistruct.BR.qvals;
qq(qq == 0) = 1e-20;
guistruct.BR.qvals(qq(qq == 0)) = 1e-20;
qval = -log10(qq);

% PCA
%scor = guistruct.scores;
%load = guistruct.loadings(:,1);

% LDA
scor = guistruct.scores;
wght = guistruct.weights(:,1);

figure; 

ax(1) = subplot(1,3,1); hold on;
[~,~,i] = unique(op.histID(idx,:));
cols = jet(numel(gnam));
for n = 1:numel(gnam)
    fx = gidx == n;
    h(n) = scatter(scor(fx,1),scor(fx,2),80,cols(n,:),'o','filled');
end
legend(h,gnam,'Location','North');
xlabel('LDA CV1','FontSize',16,'FontWeight','bold');
ylabel('LDA CV2','FontSize',16,'FontWeight','bold');
set(gca,'FontSize',14);

ax(2) = subplot(1,3,2);
scatter(wght,qval,30,'blue','o','filled');
xlabel('LDA Weight 1','FontSize',16,'FontWeight','bold');
ylabel('-log_{10} (q-val)','FontSize',16,'FontWeight','bold');
set(gca,'FontSize',14);
ylim([-0.1 20.1]);


ax(3) = subplot(1,3,3);
boxplot(guistruct.X(:,1),gidx,'Labels',gnam);

    
set(ax(2),'ButtonDownFcn',{@ccb,wght,qval,...
    guistruct.ppm,...
    guistruct.X,...
    ax(3),...
    gnam,...
    gidx});
    




end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ccb(src,event,x,y,dataX,dataY,ax,gnam,gidx)

coor = get(src,'CurrentPoint');

xc = coor(1,1);
yc = coor(1,2);

xd = (x - xc) .^ 2;
yd = (y - yc) .^ 2;

xd = xd / max(xd);
yd = yd / max(yd);

[td,id] = min(xd + yd');

%[x(id) y(id)];

% Now let's boxplot the data...
axes(ax);
cla reset;
boxplot(dataY(:,id),gidx,'Labels',gnam);
title(['m/z = ' sprintf('%0.4f',dataX(id)) ...
    ', q = ' sprintf('%0.4E',10^-(y(id)))],...
    'FontSize',16,'FontWeight','bold');
ylabel('log_{10} intensity','FontSize',16,'FontWeight','bold');
set(gca,'FontSize',14);



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
