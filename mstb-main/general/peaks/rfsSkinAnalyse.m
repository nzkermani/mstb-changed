function [ data ] = rfsSkinAnalyse(sp,data)
%rfsSkinAnalyse - perform PCA/UV of peaks in the data

% File information
fold = '/Volumes/JSM/DB/Renata/';
file = {'Skin-2C-Full.mat';'Skin-2C-NR5.mat';'Skin-2C-NR15.mat'};
numF = numel(file);

% First we just read in the data, extract annotatated parts
if isempty(data)
    [data] = importData(fold,file,numF);
    return
end

% Compare PCA
%pcaCompare(data,sp);

% Plot to show matching / missing peaks...
matchMiss(data,sp);


% Perform PCA
%pcaSmall(data(3),sp{3,3})

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = importData(fold,file,numF)
% Read in the data...

% Create structure
data = struct('file',[],'mz',[],'sp',[],'histID',[]);

for n = 1:numF
    
    % Load each file
    tmp = open([fold file{n}]);
    
    % Extract all annotated regions
    [mask2,histID,~] = desiAnnotationExtract(tmp.dpn);
    
    % Remove background
    mask2(mask2 == 15) = 0;
    
    fx = mask2 > 0;
        
    sp = tmp.dpn.d1.sp;
    sz = size(sp);
    sp = reshape(sp,[sz(1)*sz(2) sz(3)]);
    sp = sp(fx,:);
    
    fy = nanmean(sp,1) > 0;

    data(n).file = file{n};
    data(n).mz = tmp.dpn.d1.mz(fy);
    data(n).sp = sp(:,fy);
    data(n).histID = histID(fx,:);
    
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = pcaSmall(data,exc)
% Perform PCA of the smallest dataset, with all peaks, and without those
% that were not matched to the larger datasets

tmp = data.sp;
tmp = bsxfun(@rdivide,tmp,nansum(tmp,2)) * 1000;
os = nanmedian(tmp(tmp > 0));
tmp = log(tmp + os);
[l1,ss,ee] = pca(tmp,'NumComponents',2);
ee = 100 * ee / sum(ee);
titText = ['NR15 | PC1 (' sprintf('%0.1f',ee(1)) '%) v PC2 (' sprintf('%0.1f',ee(2)) '%)'];
scatterPlotNice(ss(:,1:2),data.histID,[1 0 0; 0 0 1],titText);
box on;

% Now exclude the others...
tmp = data.sp(:,~exc(:,1));
tmp = bsxfun(@rdivide,tmp,nansum(tmp,2)) * 1000;
os = nanmedian(tmp(tmp > 0));
tmp = log(tmp + os);
[l2,ss,ee] = pca(tmp,'NumComponents',2);
ee = 100 * ee / sum(ee);
titText = ['NR15 not in Full | PC1 (' sprintf('%0.1f',ee(1)) '%) v PC2 (' sprintf('%0.1f',ee(2)) '%)'];
scatterPlotNice(ss(:,1:2),data.histID,[1 0 0; 0 0 1],titText);
box on;

% Now exclude the others...
tmp = data.sp(:,~exc(:,2));
tmp = bsxfun(@rdivide,tmp,nansum(tmp,2)) * 1000;
os = nanmedian(tmp(tmp > 0));
tmp = log(tmp + os);
[l3,ss,ee] = pca(tmp,'NumComponents',2);
ee = 100 * ee / sum(ee);
titText = ['NR15 not in NR5 | PC1 (' sprintf('%0.1f',ee(1)) '%) v PC2 (' sprintf('%0.1f',ee(2)) '%)'];
scatterPlotNice(ss(:,1:2),data.histID,[1 0 0; 0 0 1],titText);
box on;

%-----------------------

figure; hold on;

stem(data.mz,l1(:,1),'k','LineWidth',2,'MarkerSize',0.01);
stem(data.mz(~exc(:,1)),l2(:,1),'r','LineWidth',2,'MarkerSize',0.01);
stem(data.mz(~exc(:,2)),l3(:,1),'b','LineWidth',2,'MarkerSize',0.01);



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function matchMiss(data,sp)

figure;
numSP = 3;
ax = zeros(numSP,1);


i = 1;
ax(i,1) = subplot(numSP,1,i); hold on;
stem(data(3).mz,nanmean(data(3).sp,1),'Color','k',...
    'LineWidth',2,'MarkerSize',0.01);
set(gca,'FontSize',14);
title('NR15 | All Peaks','FontSize',16);
box on;

% NR5 peaks not in Full
i = 2;
ax(i,1) = subplot(numSP,1,i); hold on;
inc = sp{3,3}(:,1);
t1 = data(3).sp(:,inc);
t2 = data(3).sp(:,~inc);
stem(data(3).mz(~inc),nanmean(t2,1),'Color','r',...
    'LineWidth',2,'MarkerSize',0.01);
stem(data(1).mz,-nanmean(data(1).sp,1),'Color','b',...
    'LineWidth',2,'MarkerSize',0.01);
scatter(data(3).mz(~inc),zeros(size(t2,2),1),60,'g','o','filled');
ylim([-max(nanmean(t2,1)) max(nanmean(t2,1))]);
set(gca,'FontSize',14);
title('Peaks in NR15 not Matched to Peaks in Full','FontSize',16);
box on;

% NR5 peaks not in NR15
i = 3;
ax(i,1) = subplot(numSP,1,i); hold on;
inc = sp{3,3}(:,2);
t1 = data(3).sp(:,inc);
t2 = data(3).sp(:,~inc);
stem(data(3).mz(~inc),nanmean(t2,1),'Color','r',...
    'LineWidth',2,'MarkerSize',0.01);
stem(data(2).mz,-nanmean(data(2).sp,1),'Color','b',...
    'LineWidth',2,'MarkerSize',0.01);
scatter(data(3).mz(~inc),zeros(size(t2,2),1),60,'g','o','filled');
ylim([-max(nanmean(t2,1)) max(nanmean(t2,1))]);
set(gca,'FontSize',14);
title('Peaks in NR15 not Matched to Peaks in NR5','FontSize',16);
box on;

xlabel('m/z','FontSize',16);

linkaxes(ax,'x');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pcaCompare(data,sp)
% PCA plots for peaks...


% Peaks in Full and NR15 not found in NR5 - aim to see if the noisy peaks
% that were removed are useful in any way...
fx = isnan(sp{1,3});
mz1 = data(1).mz(fx);
tmp1 = data(1).sp(:,fx);
tmp1 = bsxfun(@rdivide,tmp1,nansum(tmp1,2)) * 1000;
os = nanmedian(tmp1(tmp1 > 0));
[pqFull] = univariate(tmp1,data(1).histID,'Test','anova');
tmp1 = log(tmp1 + os);

[lf,sf,ef] = pca(tmp1,'NumComponents',2);
ef = 100 * ef / sum(ef);

% Now do peaks in NR15 not found in NR5
fx = isnan(sp{2,3});
mz2 = data(2).mz(fx);
tmp2 = data(2).sp(:,fx);
tmp2 = bsxfun(@rdivide,tmp2,nansum(tmp2,2)) * 1000;
os = nanmedian(tmp2(tmp2 > 0));
[pq15] = univariate(tmp2,data(2).histID,'Test','anova');
tmp2 = log(tmp2 + os);

[l15,s15,e15] = pca(tmp2,'NumComponents',2);
e15 = 100 * e15 / sum(e15);

% New Loadings to account for anova
ldFull = [lf(:,1) lf(:,1)];
ldFull(~(pqFull(:,2) < 0.01),2) = NaN;
ldNR15 = [l15(:,1) l15(:,1)];
ldNR15(~(pq15(:,2) < 0.01),2) = NaN;

figure;
ax(1) = subplot(2,2,1);
ax(2) = subplot(2,2,2);
ax(3) = subplot(2,2,3);
ax(4) = subplot(2,2,4);

cols = [26 150 65;215 25 28] / 255;

titText = 'Peaks from Full not in NR15';
scatterPlotNice(sf(:,1:2),data(1).histID,cols,titText,ax(1));

titText = 'Peaks from NR5 not in NR15';
scatterPlotNice(s15(:,1:2),data(2).histID,cols,titText,ax(3));

numV1 = sum(~isnan(ldFull(:,2)));
titText = ['PC1 Loadings | ' sprintf('%d',numV1) '/' sprintf('%d',size(ldFull,1)) ' variables with ANOVA q < 0.01'];    
stemPlotNice(mz1,ldFull,titText,ax(2));

numV2 = sum(~isnan(ldNR15(:,2)));
titText = ['PC1 Loadings | ' sprintf('%d',numV2) '/' sprintf('%d',size(ldNR15,1)) ' variables with ANOVA q < 0.01'];    
stemPlotNice(mz2,ldNR15,titText,ax(4));


linkaxes([ax(2) ax(4)],'x');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%