function [ output_args ] = uiMSreimg(g,origPath)
%uiMSreimg - make prettier plots than are available from the stats toolbox

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

% Determine group labels...
info.num = g.groupdata';
info.ids = g.groupIds;
%info.ids = {'Borderline','Cancer','Healthy (OV)'};
%info.ids = {'Cancer','Healthy (FT)','Healthy (OV)'};
%info.ids = {'Clear Cell','Endometrioid','Healthy (OV)','Serous'};
%info.ids = {'Healthy (OV)','Serous','Stroma'};
%info.ids = {'Healthy (FT)','Healthy (OV)','Serous','Stroma'};
%info.ids = {'Cancer','Healthy'};


% Reformat for box plot compatibility
info.all = cell(numel(info.num),1);
for n = 1:numel(info.ids)
    fx = info.num == n;
    info.all(fx) = info.ids(n);
end

% Check colours...
try
    if numel(info.ids) > size(cols,1)
        cols = g.spcolors;
    end
catch
    cols = g.spcolors;
end

% PCA information
pca.cm = [1 2];
pca.sc = g.pca.T(:,pca.cm);
pca.ld = g.pca.P(:,pca.cm);
pca.eg = g.pca.e(pca.cm);

% MMC & CV information
mmc.cm = [1 2];
mmc.sc = g.cv.T(:,mmc.cm);
mmc.ld = g.weights(:,mmc.cm);%g.cv.weights(:,mmc.cm);
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



% Need to do the confusion matrix...
makeConfMat(mmc.conf.acc,mmc.conf.smpls,cols,order);
graphFormat([path 'ConfMat'],'png');

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

% Loadings/weights plot
l1 = 'm/z';
l2 = 'Component Loading';
loadingsPlot(g.ppm,pca.ld(:,1),l1,l2);
savefig([path 'PC1-Loadings.fig']);
graphFormat([path 'PC1-Loadings'],'png');
loadingsPlot(g.ppm,pca.ld(:,2),l1,l2);
savefig([path 'PC2-Loadings.fig']);
graphFormat([path 'PC2-Loadings'],'png');

% Loadings/weights plot
l1 = 'm/z';
l2 = 'Component Loading';
loadingsPlot(g.ppm,mmc.ld(:,1),l1,l2);
savefig([path 'MMC1-Loadings.fig']);
graphFormat([path 'MMC1-Loadings'],'png');
loadingsPlot(g.ppm,mmc.ld(:,2),l1,l2);
savefig([path 'MMC2-Loadings.fig']);
graphFormat([path 'MMC2-Loadings'],'png');

% Now do the PCA scatter plot
l1 = ['PC' int2str(pca.cm(1)) ', ' sprintf('%d',round(pca.eg(1))) '%  '];
l2 = ['PC' int2str(pca.cm(2)) ', ' sprintf('%d',round(pca.eg(2))) '%  '];
scatterPlot(pca.sc(:,1),pca.sc(:,2),info.num,info.ids,cols,l1,l2,order);
graphFormat([path 'PCA-Scores'],'png');

% Now do the MMC scatter plot
l1 = ['LV' int2str(1)];
l2 = ['LV' int2str(2)];
scatterPlot(mmc.sc(:,1),mmc.sc(:,2),info.num,info.ids,cols,l1,l2,order);
graphFormat([path 'MMC-Scores'],'png');

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


xlabel(xl,'FontSize',18,'FontWeight','bold');
ylabel(yl,'FontSize',18,'FontWeight','bold');


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadingsPlot(x,y,l1,l2)

figure('Position',[700 100 700 350]);
hold on;

stem(x,y,'k','MarkerSize',0.1);

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
function makeConfMat(prc,qty,cols,order)

% Size of the matrices
sz = size(prc,1);

if isempty(order)
    order = 1:sz
end


% We need to try to sort the order of the confusion matrices
prc2 = zeros(size(prc));
qty2 = zeros(size(qty));
for n = 1:sz
    for r = 1:sz
        prc2(n,r) = prc(order(n),order(r));
        qty2(n,r) = qty(order(n),order(r));
    end
end
prc = prc2;
qty = qty2;

cmap = flipud(gray(100));
fs = 18;

figure; hold on;

imagesc(prc);
colormap(cmap);
caxis([0 100]);

set(gca,'YDir','reverse','XTick',[],'YTick',[]);

scatter(1:sz,zeros(sz,1)+0.5,200,cols(order,:),'o','filled',...
    'MarkerEdgeColor',[0.5 0.5 0.5]);

scatter(zeros(sz,1)+0.5,1:sz,200,cols(order,:),'o','filled',...
    'MarkerEdgeColor',[0.5 0.5 0.5]);

for n = 1:sz
    
    for r = 1:sz
        
        % Skip empty pixels
        if qty(r,n) == 0
            continue;
        end
        
        
        if prc(r,n) > 60
            tc = 'white';
        else
            tc = 'black';
        end
        
        tl = [sprintf('%0.1f',prc(r,n)) '%' char(10) sprintf('%d', qty(r,n))];
        
        text(n,r,tl,...
            'Color',tc,...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...
            'FontSize',fs,...
            'FontWeight','bold');
        
    end
    
end


ylabel('Actual Class   ','FontSize',20,'FontWeight','bold');
title('Predicted Class   ','FontSize',20,'FontWeight','bold');
axis tight

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%