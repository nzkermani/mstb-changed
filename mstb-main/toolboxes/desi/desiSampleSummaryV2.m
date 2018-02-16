function [fig] = desiSampleSummaryV2(dpn,mzVals)
%desiSampleSummary

% Various parameters
threshold = 0.985;
mzRange = [600 1000];
ppmTol = 5;


% Determine colours for the various annotated regions
[cols,unqH] = getAnnotationColours(dpn.anno);

% Indices of the training set
[dx1,~,anno] = dpnPrepFullMS(dpn,mzRange,[]);

% Find the variable mz values
numM = numel(dx1.mz);
mzMask = zeros(numel(mzVals),numM);
ppmDev = ppmTol * mzVals / 1e6;
for n = 1:numel(mzVals)

    mmm = mzVals(n) + [-ppmDev(n) ppmDev(n)];
    mzMask(n,:) = dx1.mz > mmm(1) & dx1.mz < mmm(2);    

end
mzMask = sum(mzMask,1) > 0;

% Normalise / transform the data
dx1.sp = bsxfun(@rdivide,dx1.sp,nansum(dx1.sp,2));
dx1.sp = 1e4 * dx1.sp / nanmax(dx1.sp(:));
dx1.sp(isnan(dx1.sp)) = 0;

% Now log transform
os = nanmedian(dx1.sp(dx1.sp > 0));
dx1.sp = log(dx1.sp + os);

% Run the MMC
[unam,anno,predID] = runMMC(dx1,anno,threshold);

% Reshape the unambiguous predictions into an image then coloured according
% to the annotation colours. Then plot it...
sz = size(dpn.d1.sp);
clImg = reshape(unam,[sz(1) sz(2)]);
clImg = kmeans2rgb(clImg,cols);




% Draw a figure
fname = dpn.file.nam;
[fig] = drawFigure(fname);
drawImages(dpn,clImg,fig,unqH,cols);

% Box plots for all predicted pixels
fx = unam > 0;
jsmBoxPlotMulti(dx1.sp(fx,mzMask),predID(fx),dx1.mz(mzMask),...
    'Cols',cols,...
    'Orientation',fig.bpPred);
set(fig.bpPred,'FontSize',12);

% Box plots for only the annotated pixels
fy = anno.mask2 > 0;
jsmBoxPlotMulti(dx1.sp(fy,mzMask),anno.histID(fy),dx1.mz(mzMask),...
    'Cols',cols,...
    'Orientation',fig.bpAnno);
set(fig.bpAnno,'FontSize',12);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = drawFigure(fname)

% Draw a figure
fig.fig = figure('Units','normalized','Position',[0.1 0.5 0.35 0.88],...
    'Color','white');

% Title...
%fname = dpn.file.nam;
dot = strfind(fname,'.');
fname = fname(1:dot(end)-1);
uicontrol('Parent',fig.fig,'Units','normalized',...
    'Position',[0.03 0.97 0.5 0.02],...
    'Style','text',...
    'String',fname,...
    'FontSize',20,...
    'FontWeight','bold',...
    'HorizontalAlignment','left',...
    'BackgroundColor','white');

% Axes for optical image and MS images with annotations...
fig.axOpt = axes('Parent',fig.fig,'Units','normalized','Position',[0.05 0.7 0.4 0.25]);
fig.axMSI = axes('Parent',fig.fig,'Units','normalized','Position',[0.55 0.7 0.4 0.25]);

fig.axGrp = axes('Parent',fig.fig,'Units','normalized','Position',[0.05 0.40 0.4 0.25]);
fig.axPrd = axes('Parent',fig.fig,'Units','normalized','Position',[0.55 0.40 0.4 0.25]);

fig.bpAnno = axes('Parent',fig.fig,'Units','normalized','Position',[0.05 0.1 0.4 0.25]);
fig.bpPred = axes('Parent',fig.fig,'Units','normalized','Position',[0.55 0.1 0.4 0.25]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [unam,anno,predID] = runMMC(dx1,anno,threshold)
% Run MMC and determine the pixels that are unambiguously annotated

% Indices of annotated pixels
idx = anno.mask2 > 0;

% What is the mean spectrum?
mean1 = nanmean(dx1.sp(idx,:),1);

% Subtract the mean.
tmp1 = bsxfun(@minus,dx1.sp(idx,:),mean1);

% Run the OAA-MMC function
[B1,dx1.load] = oaaMMC(tmp1,anno.histID(idx),2);

% Calculate probabilities for pixels
[~,Pr1,~] = calcScoresSpectra(dx1.sp,mean1,dx1.load,B1);

% Save probabilities
dx1.data = Pr1;

% Find unambiguously classified pixels, and return a matrix of class
% membership
unam = Pr1 > threshold;
su = sum(unam,2);
unam(su > 1,:) = 0;
unam = bsxfun(@times,unam,1:size(unam,2));
unam = max(unam,[],2);

% Rather than numbers, return IDs in text string format
predID = cell(size(anno.histID));
[unqH,~,~] = unique(anno.histID(anno.mask2 > 0));
numH = numel(unqH);
for n = 1:numH
    fx = unam == n;
    predID(fx,1) = unqH(n);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cols,unqH] = getAnnotationColours(anno)
% Determine the colours for the annotated regions

numA = size(anno,1);
[unqH,~,~] = unique(anno(:,5));
numH = numel(unqH);
cols = zeros(numH,3);
for n = 1:numA
    
    % Save the colour information
    fx = strcmp(unqH,anno{n,5});
    cols(fx,:) = anno{n,3};    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawImages(dpn,clImg,fig,unqH,cols)
% Draw the two images in the right places, with annotation patches

% Optical image
imagesc(fig.axOpt,dpn.opt.coreg);
set(fig.axOpt,'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(dpn.opt.coreg,2)],...
    'YLim',[1 size(dpn.opt.coreg,1)],...
    'LineWidth',5,...
    'TickLength',[0 0]);

% Add in the annotations...
numA = size(dpn.anno,1);
axes(fig.axOpt); hold on;
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
imagesc(fig.axMSI,nansum(dpn.d1.sp,3));
set(fig.axMSI,'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(dpn.d1.sp,2)],...
    'YLim',[1 size(dpn.d1.sp,1)],...
    'LineWidth',5,...
    'TickLength',[0 0]);

% Add in the annotations...
numA = size(dpn.anno,1);
axes(fig.axMSI); hold on;
for n = 1:numA
    redrawPatch(dpn.anno{n,8},dpn.anno{n,9},dpn.anno{n,3});    
    
    % Save the colour information
    %fx = strcmp(unqH,dpn.anno{n,5});
    %cols(fx,:) = dpn.anno{n,3};    
end

% Now do the annotations for show
axes(fig.axGrp); hold on;
for n = 1:numel(unqH)
    
    scatter(0,-n,300,cols(n,:),'o','filled','MarkerEdgeColor','k');
    
    text(0.2,-n,unqH{n},'FontSize',20,'FontWeight','bold',...
        'HorizontalAlignment','left');
    
end
xlim(fig.axGrp,[-0.05 2]);
ylim(fig.axGrp,[(-numel(unqH))-0.5 -0.5]);
axis off;

% Predicted image...
imagesc(fig.axPrd,clImg);
set(fig.axPrd,'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(clImg,2)],...
    'YLim',[1 size(clImg,1)],...
    'LineWidth',5,...
    'TickLength',[0 0]);

end
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