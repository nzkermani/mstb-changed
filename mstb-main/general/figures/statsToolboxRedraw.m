function statsToolboxRedraw(g,origPath)
%uiMSreimg - make prettier plots than are available from the stats toolbox

% Which plots to draw? Specify these here
opts.doPCA = false;
opts.doMMC = false;
opts.doCM  = false;
opts.avgsp = false;
opts.univar= true;

% Define some additional parameters if necessary?
opts.order = [];
opts.pcs   = [1 2];
opts.mmcs  = [1 2];
opts.pval  = 0.01;

% Want to define the colours yourself? Leave [] otherwise
opts.cols  = [];
%opts.cols = [227 23 168; 108 245 149] / 256; % cut
%opts.cols = [227 94 23; 12 156 115] / 256; %& coag
%opts.cols = [237 29 19; 12 133 34] / 256; % cut/coag comb
%opts.cols = [227 94 23; 227 23 168; 12 156 115; 108 245 149] / 256; % cut/coag sep
opts.cols = [237 29 19; 0 153 255] / 256; % fibro v canc
%
%
% Change only the stuff above
%
%

% Determine group labels...
info.num = g.groupdata';
info.ids = g.groupIds;

% Reformat for box plot compatibility
info.all = cell(numel(info.num),1);
for n = 1:numel(info.ids)
    fx = info.num == n;
    info.all(fx) = info.ids(n);
end

% Check colours...
try
    if numel(info.ids) > size(opts.cols,1)
        opts.cols = g.spcolors;
    end
catch
    opts.cols = g.spcolors;
end

% Extract the PCA information from the guidata structure
pca.sc = g.pca.T(:,opts.pcs);
pca.ld = g.pca.P(:,opts.pcs);
pca.eg = g.pca.e(opts.pcs);

% MMC & CV information
mmc.sc = g.cv.T(:,opts.mmcs);
mmc.ld = g.weights(:,opts.mmcs);%g.cv.weights(:,opts.mmcs);
mmc.conf = g.cv.confMat;
mmc.ac = g.cv.accuracy(g.cv.nComps);

% Here we need a file name / folder generator...
if nargin == 2
    tmp = strfind(origPath,'.');
    path = [origPath(1:tmp(end)-1) '/'];    
else
    path = [pwd filesep 'Images' filesep];    
end

if ~exist(path,'dir')
        mkdir(path);
end


% Confusion matrix   
if opts.doCM
    
    % Plot and save
    makeConfMat(mmc.conf.acc,mmc.conf.smpls,opts.cols,opts.order);
    graphFormat([path 'ConfMat'],'png');
end

% Now do the PCA scatter plot
if opts.doPCA
    
    % Labels
    l1 = ['PC' int2str(opts.pcs(1)) ', ' sprintf('%d',round(pca.eg(1))) '%  '];
    l2 = ['PC' int2str(opts.pcs(2)) ', ' sprintf('%d',round(pca.eg(2))) '%  '];
    
    % Plot and save
    scatterPlot(pca.sc(:,1),pca.sc(:,2),...
        info.num,info.ids,opts.cols,l1,l2,opts.order);
    
    graphFormat([path 'PCA-Scores'],'png');
end

% Now do the MMC scatter plot
if opts.doMMC
        
    % Labels
    l1 = ['LV' int2str(opts.mmcs(1))];
    
    % Only show a single LV if 2 classes
    if size(g.groupIds,1) > 2
        l2 = ['LV' int2str(opts.mmcs(2))];
    else
        mmc.sc(:,2) = randn(size(mmc.sc,1),1);    
        l2 = 'Random Index';
    end
    
    % Plot and save
    scatterPlot(mmc.sc(:,1),mmc.sc(:,2),...
        info.num,info.ids,opts.cols,l1,l2,opts.order);
    
    graphFormat([path 'MMC-Scores'],'png');
end

% What about average spectra for these groups?
if opts.avgsp
    
    plotAverageSpectra(g.ppm,g.Sp,info.num,info.ids,opts.cols)
    
    graphFormat([path 'Average-Spectra'],'png');
end    


if opts.univar
    
    plotUnivariate(g.ppm,g.Sp,info.num,info.ids,opts.cols,opts.pval,path);

end




return

order = [];

% Definition of colours...
%cols = [232 217 51; 255 64 0; 64 168 35]/256; %borderline,cancer,healthyOV
%order = [3 1 2];

%cols = [191 65 23; 19 179 232; 64 168 35]/256; %cancer,ft,ov
%order = [3 2 1];

%cols = [227 64 137; 255 169 20; 64 168 35; 149 23 191]/256; %clear,endo,norm,serous
%order = [3 1 2 4];

%cols = [64 168 35; 149 23 191; 232 217 51]/256; % ov,ser,stroma
%order = [1 3 2];

%cols = [19 179 232; 64 168 35; 149 23 191; 232 217 51]/256; % ft,ov,ser,stroma
%order = [1 2 4 3];


%cols = [255 64 0; 64 168 35]/256; %cancer,healthyOV
%order = [2 1];

%info.ids = {'Borderline','Cancer','Healthy (OV)'};
%info.ids = {'Cancer','Healthy (FT)','Healthy (OV)'};
%info.ids = {'Clear Cell','Endometrioid','Healthy (OV)','Serous'};
%info.ids = {'Healthy (OV)','Serous','Stroma'};
%info.ids = {'Healthy (FT)','Healthy (OV)','Serous','Stroma'};
%info.ids = {'Cancer','Healthy'};









% % Boxplots...
% [~,idx] = max(pca.ld(:,1));
% jsmBoxPlot(g.Sp(:,idx),info.all,cols,order);
% fn = ['BP-mz-' sprintf('%0.3f',g.ppm(idx))];
% tmp = strfind(fn,'.');
% fn(tmp) = '-';
% graphFormat([path fn],'png');
% 
% [~,idx] = min(pca.ld(:,1));
% jsmBoxPlot(g.Sp(:,idx),info.all,cols,order);
% fn = ['BP-mz-' sprintf('%0.3f',g.ppm(idx))];
% tmp = strfind(fn,'.');
% fn(tmp) = '-';
% graphFormat([path fn],'png');
% 
% [~,idx] = max(mmc.ld(:,1));
% jsmBoxPlot(g.Sp(:,idx),info.all,cols,order);
% fn = ['BP-mz-' sprintf('%0.3f',g.ppm(idx))];
% tmp = strfind(fn,'.');
% fn(tmp) = '-';
% graphFormat([path fn],'png');
% 
% [~,idx] = min(mmc.ld(:,1));
% jsmBoxPlot(g.Sp(:,idx),info.all,cols,order);
% fn = ['BP-mz-' sprintf('%0.3f',g.ppm(idx))];
% tmp = strfind(fn,'.');
% fn(tmp) = '-';
% graphFormat([path fn],'png');

% % Loadings/weights plot
% l1 = 'm/z';
% l2 = 'Component Loading';
% loadingsPlot(g.ppm,pca.ld(:,1),l1,l2);
% savefig([path 'PC1-Loadings.fig']);
% graphFormat([path 'PC1-Loadings'],'png');
% loadingsPlot(g.ppm,pca.ld(:,2),l1,l2);
% savefig([path 'PC2-Loadings.fig']);
% graphFormat([path 'PC2-Loadings'],'png');
% 
% % Loadings/weights plot
% l1 = 'm/z';
% l2 = 'Component Loading';
% loadingsPlot(g.ppm,mmc.ld(:,1),l1,l2);
% savefig([path 'MMC1-Loadings.fig']);
% graphFormat([path 'MMC1-Loadings'],'png');
% loadingsPlot(g.ppm,mmc.ld(:,2),l1,l2);
% savefig([path 'MMC2-Loadings.fig']);
% graphFormat([path 'MMC2-Loadings'],'png');


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scatterPlot(x,y,grp,names,cols,xl,yl,ord)

ms = 120;

numG = size(cols,1);

if nargin < 8
    ord = 1:numG;
end

if isempty(ord)
    ord = 1:numG;
end

figure('Position',[100 100 800 600]); hold on;

h = zeros(numG,1);
for n = 1:numG
    
    fx = grp == ord(n);
    
    h(n) = scatter(x(fx),y(fx),ms,cols(ord(n),:),'o','filled',...
        'MarkerEdgeColor',[0.5 0.5 0.5]);
    
end

% Overall font size
set(gca,'FontSize',16);

% Legend placement and location
leg = legend(h,names(ord),'Location','NorthEastOutside');
m = findobj(leg,'Type','patch');
set(m,'MarkerSize',sqrt(ms));

box on;

% Axes limits
xL = [min(x) max(x)];
yL = [min(y) max(y)];
xS = (xL(2) - xL(1)) * 0.025;
yS = (yL(2) - yL(1)) * 0.025;
xlim([xL(1)-xS xL(2)+xS]);
ylim([yL(1)-yS yL(2)+yS]);

% Labels and formatting
xlabel(xl,'FontSize',18,'FontWeight','bold');
ylabel(yl,'FontSize',18,'FontWeight','bold');
axis square

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadingsPlot(x,y,l1,l2)

figure('Position',[700 100 700 350]);
hold on;

stem(x,y,'k','MarkerSize',0);

% Determine better axes limits
yl = [min(y) max(y)];
ylim([yl(1)-0.01 yl(2)+0.01]);

% Overall font size
set(gca,'FontSize',16);

box on;

% Axes labels
xlabel(l1,'FontSize',18,'FontWeight','bold');
ylabel(l2,'FontSize',18,'FontWeight','bold');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotAverageSpectra(mz,sp,grp,names,cols)

% How many groups?
numG = size(cols,1);

% Draw a figure
figure('Position',[100 100 1200 500]); hold on;

% Figure handles
h = zeros(numG,1);

% Loop through each group
for n = 1:numG
    
    fx = grp == n;
    
    % Calculate the mean spectrum
    mns = nanmean(sp(fx,:),1);
    [m1,m2] = insertZeros(mz,mns);
    
    % If numG == 2 and this is two, then reverse it
    if numG == 2 && n == 2
        m2 = m2 * -1;
    elseif numG == 4 && n >= 3
        m2 = m2 * -1;
    end
    
    h(n) = plot(m1,m2,'Color',cols(n,:),'LineWidth',1);
    
end

% Overall font size
set(gca,'FontSize',16);

% Legend placement and location
leg = legend(h,names,'Location','NorthEastOutside');
%m = findobj(leg,'Type','patch');
%set(m,'MarkerSize',sqrt(ms));

box on;

% Axes limits
% xL = [min(x) max(x)];
% yL = [min(y) max(y)];
% xS = (xL(2) - xL(1)) * 0.025;
% yS = (yL(2) - yL(1)) * 0.025;
xlim([600 1000]);
% ylim([yL(1)-yS yL(2)+yS]);

% Labels and formatting
xlabel('m/z   ','FontSize',18,'FontWeight','bold','FontAngle','italic');
%ylabel('Intensity','FontSize',18,'FontWeight','bold');
set(gca,'YTick',[]);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotUnivariate(mz,sp,grp,names,cols,pThresh,path)
% Perform univariate statistics and then plot the most significant

numV = numel(mz);
numG = numel(grp);

pvals = zeros(numV,3);
pvals(:,1) = 1:numV;

for n = 1:numV
    
    % Do Kruskal Wallis
    [pvals(n,2)] = kruskalwallis(sp(:,n),grp,'off');
    
end

% Q-value correction
tmp = getBHYqVls(pvals(:,2)',0.001);
pvals(:,3) = tmp';

% Sort pqvals
pvals = sortrows(pvals,3);

idx = pvals(:,3) < pThresh;

% Now make box plots for the x lowest q values
for n = 1:min([30 sum(idx)])
    
    % Variable index
    i = pvals(n,1);
    
    % What about a title?
    tit = ['m/z = ' sprintf('%0.2f',mz(i))];
    
    jsmBoxPlot(sp(:,i),grp,...
        'Labels',names,...
        'Colours',cols,...
        'Legend',false,...
        'Orientation','vertical',...
        'Title',tit);
    
    
    graphFormat([path 'Boxplot-' int2str(i) '-mz' int2str(mz(i))],'png');
    close
    

end
    
    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

