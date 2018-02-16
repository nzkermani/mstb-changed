function dbStats(data)
% dbStats - quick function for analysing the results from multiple sample
% analysis.

% OPTIONS
opts.mzRange    = [100 1000];
opts.normalise  = 'tic';
opts.minSamples = 0;
opts.log        = false;


% Trim variables
[data] = mzTrim(data,opts.mzRange,opts.minSamples);

% Normalise
[data] = normalise(data,opts.normalise);

% Spectral average plot
spectralAverage(data);

% Log data
[data,~] = logData(data,opts.log);

% PCA
msaPCA(data);

% MMC
msaMMC(data);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = mzTrim(data,mzRange,minSamp)
% Remove variables outside of the specified range, and any which don't
% appear in the minimum number of samples

% Find frequency of variables in each file
[unq,~,ind] = unique(data.meta.fileID);
numF = numel(unq);

freq = zeros(numF,numel(data.mz));
for n = 1:numF
    
    idx = ind == n;
    %freq(n,:) = sum(data.sp(idx,:) > 0,1);
    freq(n,:) = sum(data.sp(idx,:) > 0,1) >= 1;%(0.5 * sum(idx));
    
end

% Sum freq over the rows
freq = sum(freq,1) > minSamp;
mask = data.mz >= min(mzRange) & data.mz <= max(mzRange);

% Combine the masks
comb = mask & freq';

data.mz = data.mz(comb);
data.sp = data.sp(:,comb);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = normalise(data,method)

switch method
    
    case 'tic'
        tic = nansum(data.sp,2);
        data.sp = bsxfun(@rdivide,data.sp,tic) * 1e4;
        
    otherwise
        disp('No normalisation performed');
end

fx = tic == 0;
if sum(fx > 0)
    data.sp = data.sp(~fx,:);
    data.meta.histID = data.meta.histID(~fx,:);
    data.meta.fileID = data.meta.fileID(~fx,:);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,os] = logData(data,method)

if method    
    os = nanmedian(data.sp(data.sp > 0));    
    data.sp = log(data.sp + os);
else
    os = [];
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spectralAverage(data)
% Simple average spectrum of all groups

new = classMany2One([data.meta.histID data.meta.fileID]);
groupPlot(data.mz,data.sp,new,'average');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msaPCA(data)
% PCA

[pc.ll,pc.ss,pc.ee] = pca(data.sp,'NumComponents',10);
pc.ee = 100 * pc.ee / sum(pc.ee);

% Axes labels
labX = ['PC1 (' sprintf('%0.1f',pc.ee(1)) '%)'];
labY = ['PC2 (' sprintf('%0.1f',pc.ee(2)) '%)'];

% Plots
[~,~] = scatterPlot(pc.ss(:,1:2),data.meta.histID,[],labX,labY);
[~,~] = scatterPlot(pc.ss(:,1:2),data.meta.fileID,[],labX,labY);

figure; hold on;
stem(data.mz,pc.ll(:,1),'b');
stem(data.mz,pc.ll(:,2),'r');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msaMMC(data)

% Run MMC with k = 10 fold cross validation.  Note that this isn't truly
% unbiased, but it is the best for now
[mmc.ss,mmc.cm,mmc.pred,mmc.unq,mmc.ll] = kfoldMMC(data.sp,data.meta.histID,10);

[~,~] = scatterPlot(mmc.ss(:,1:2),data.meta.histID,[],'LV1','LV2');

makeConfMat([],mmc.cm,[],[],mmc.unq);

figure; hold on;
stem(data.mz,mmc.ll(:,1),'b');
stem(data.mz,mmc.ll(:,2),'r');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stuff


%% Difference spectrum
figure; hold on;
idx = strcmp(histID,'Control');
df = bsxfun(@minus,nanmean(tic(idx,:),1),nanmean(tic(~idx,:),1));
gp1 = stem(mz,df,'Color',[0 0 1],'LineWidth',2,...
    'MarkerSize',0.001);
legend(gp1,{'Control - Lesion'});
box on;
xlabel('m/z','FontSize',16,'FontWeight','bold');
ylabel('Intensity','FontSize',16,'FontWeight','bold');
set(gca,'FontSize',14);
pubFig([tag 'Group-Difference'],'eps',200);
xlim([600 max(mz)]);
pubFig([tag 'Group-Difference-600-1000'],'eps',200);


%% Univariate statistics
% Use univariate statistics to identify features that are discriminatory
% between lesion and control regions (i.e. ID1 and ID5)
[pq] = univariate(tic,histID,'Test','anova');
pq(isnan(pq)) = 1;

%% MMC classification
[mmc.ss,mmc.cm,mmc.pred,mmc.unq,mmc.ll] = kfoldMMC(tic,histID,10);
%[mmd.ss,mmd.cm,mmd.pred,mmd.unq,mmd.ll] = kfoldMMC(tic,data.meta.fileID,10);

%% MMC plots
[~,~] = scatterPlot(mmc.ss(:,1:2),histID,[],'LV1','LV2');
pubFig([tag 'MMC-Scores'],'eps',200);

makeConfMat([],mmc.cm,[],[],mmc.unq);
pubFig([tag 'MMC-ConfMat'],'eps',200);

qvc = -log10(pq(:,2));
qvc(isinf(qvc)) = max(qvc);
[cfig,cax] = plotXYC(mz,mmc.ll(:,1),qvc);
box on;
xlabel('m/z','FontSize',16,'FontWeight','bold');
ylabel('LV1 Weight','FontSize',16,'FontWeight','bold');
set(cax,'FontSize',14);
c = colorbar;
colormap(redgreencmap);
set(c,'FontSize',14);
ylabel(c,'-log_{10}(q)','FontSize',16,'FontWeight','bold');
pubFig([tag 'MMC-Weights'],'eps',200);
xlim([600 max(mz)]);
pubFig([tag 'MMC-Weights-600-1000'],'eps',200);



%% Determine the q values of better variables
pqMask = pq(:,2) < 0.01;
topVar = [find(pqMask) mz(pqMask)' pq(pqMask,2)];
topVar = sortrows(topVar,3);
for n = 1:size(topVar,1)
    
    jsmBoxPlot(tic(:,topVar(n,1)),histID,...
        'Orientation','vertical',...
        'Order',[1 2],...
        'Labels',mmc.unq);
    pubFig([tag 'BoxPlot-mz' sprintf('%d',topVar(n,2))],'eps',200);
    close(gcf);
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%