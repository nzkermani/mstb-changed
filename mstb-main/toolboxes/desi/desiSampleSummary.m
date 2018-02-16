function [fig] = desiSampleSummary(dpn)
%desiSampleSummary

% Draw a figure
fig = figure('Units','normalized','Position',[0.1 0.5 0.35 0.88],...
    'Color','white');

% Title...
fname = dpn.file.nam;
dot = strfind(fname,'.');
fname = fname(1:dot(end)-1);
uicontrol('Parent',fig,'Units','normalized',...
    'Position',[0.03 0.97 0.5 0.02],...
    'Style','text',...
    'String',fname,...
    'FontSize',20,...
    'FontWeight','bold',...
    'HorizontalAlignment','left',...
    'BackgroundColor','white');

% Axes for optical image and MS images with annotations...
axOpt = axes('Parent',fig,'Units','normalized','Position',[0.05 0.7 0.4 0.25]);
axMSI = axes('Parent',fig,'Units','normalized','Position',[0.55 0.7 0.4 0.25]);

axPCA = axes('Parent',fig,'Units','normalized','Position',[0.05 0.40 0.4 0.25]);
axGrp = axes('Parent',fig,'Units','normalized','Position',[0.55 0.40 0.4 0.25]);

axMMC = axes('Parent',fig,'Units','normalized','Position',[0.05 0.1 0.4 0.25]);
axCon = axes('Parent',fig,'Units','normalized','Position',[0.55 0.1 0.4 0.25]);

%return


% Extract annotated pixels
[mask,histID,~] = desiAnnotationExtract(dpn);
fx = mask > 0;

% Reshape sp matrix
sz = size(dpn.d1.sp);
sp = reshape(dpn.d1.sp,[sz(1)*sz(2) sz(3)]);

% Trim out to leave only annotated pixels
sp = sp(fx,:);
histID = histID(fx,:);
[unqH,~,~] = unique(histID);
numH = numel(unqH);
cols = zeros(numH,3);

% Normalise the data using predefined settings
sp = bsxfun(@rdivide,sp,nansum(sp,2));
sp = 1e4 * sp / nanmax(sp(:));

% Now log transform
os = nanmedian(sp(sp > 0));
sp = log(sp + os);

% Run PCA
[~,pcaSc,~] = pca(sp,'NumComponents',2);

% Run MMC with LOOCV
cvp = 1:numel(histID);
[mmcSc,cm,preds,grps,ll] = lpoMMC(sp,histID,cvp');
%[mmcSc,cm,preds,grps,~] = kfoldMMC(sp,histID,10);





% Now draw some of the stuff...

% Optical image
imagesc(axOpt,dpn.opt.coreg);
set(axOpt,'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(dpn.opt.coreg,2)],...
    'YLim',[1 size(dpn.opt.coreg,1)],...
    'LineWidth',5,...
    'TickLength',[0 0]);

% Add in the annotations...
numA = size(dpn.anno,1);
axes(axOpt); hold on;
for n = 1:numA
    x = dpn.anno{n,6};
    y = dpn.anno{n,7};
    x = [x(1) x(2) x(2) x(1)];
    y = [y(1) y(1) y(2) y(2)];
    if x(1) == x(2) && y(1) == y(3)
        x = x(1);
        y = y(1);
    end
    redrawPatch(x,y,dpn.anno{n,3});
end
    
    

% MS image...
imagesc(axMSI,nansum(dpn.d1.sp,3));
set(axMSI,'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(dpn.d1.sp,2)],...
    'YLim',[1 size(dpn.d1.sp,1)],...
    'LineWidth',5,...
    'TickLength',[0 0]);

% Add in the annotations...
numA = size(dpn.anno,1);
axes(axMSI); hold on;
for n = 1:numA
    redrawPatch(dpn.anno{n,8},dpn.anno{n,9},dpn.anno{n,3});    
    
    % Save the colour information
    fx = strcmp(unqH,dpn.anno{n,5});
    cols(fx,:) = dpn.anno{n,3};    
end

% Now do the annotations for show
axes(axGrp); hold on;
for n = 1:numH
    
    scatter(0,-n,300,cols(n,:),'o','filled','MarkerEdgeColor','k');
    
    text(0.2,-n,unqH{n},'FontSize',20,'FontWeight','bold',...
        'HorizontalAlignment','left');
    
end
xlim(axGrp,[-0.05 2]);
ylim(axGrp,[(-numH)-0.5 -0.5]);
axis off;



% PCA scores
scatterPlotNice(pcaSc(:,1:2),histID,cols,'PCA',axPCA);

% MMC, with incorrect predictions
scatterPlotNice(mmcSc(:,1:2),histID,cols,'MMC',axMMC);
fx = preds(:,1) ~= preds(:,2);
scatter(axMMC,mmcSc(fx,1),mmcSc(fx,2),250,'k','o');

% Confusion matrix
if numH < 4
    opts.mS = 300;
    opts.fS = 24;
else
    opts.mS = 250;
    opts.fS = 18;
end
opts.cols = cols;
xyz.cm.cm = cm;
xyz.cm.names = grps;
statsPlotConfMat(xyz,axCon,opts);
set(axCon,'XTick',[],...
    'YTick',[],...
    'LineWidth',5,...
    'TickLength',[0 0]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = redrawPatch(x,y,col)

if numel(x) == 1 && numel(y) == 1
    % Then this is a single scatter point of annotation
    h = scatter(x(1),y(1),200,col,'o','filled',...
        'MarkerEdgeColor','w');
    
else
    
    % Check that we have the annotations in the correct order for optical
    % and MS annotations
    x = [min(x) max(x) max(x) min(x)];
    y = [min(y) min(y) max(y) max(y)];

    % Then multiple points, so we patch them together
    x = x + [-0.5 0.5 0.5 -0.5];
    y = y + [-0.5 -0.5 0.5 0.5];
    
    h = patch(x,y,...
        col,...
        'EdgeColor',col,...
        'FaceColor',col,...
        'FaceAlpha',0.4,...
        'LineWidth',3);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%