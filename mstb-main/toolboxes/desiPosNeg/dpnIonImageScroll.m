function dpnIonImageScroll(src,~,fig,mode,man)
% Change the ion image on display. Either through a mouse scroll, probs too
% hard, or instead via a RGB channel of three ions

% Get the guidata...
dpn = guidata(fig.fig);

% Change the layout if required...
if strcmp(get(fig.tb.layout,'State'),'off')
    set(fig.tb.layout,'State','on');
    desiChangeLayout(fig.tb.layout,[],fig);
end

% Determine if we should log offset or suppress
doLog = get(man.doLog,'Value');
doSup = get(man.suppress,'Value');

% First get the 3 ions that were selected
numV = numel(dpn.d1.mz);
mask = zeros(numV,1);
for n = 1:3
    
    % Get the selected ions
    fx = get(man.list(n),'Value') - 1;
    
    % Remove any zero values as these mean leave the thing blank
    fx(fx == 0) = [];   
    
    % Add ions to the list to be summed later on
    mask(fx,n) = 1;
end

% Determine the image to be plotted...
sz = size(dpn.d1.sp);
new = zeros(sz(1),sz(2),3);
for n = 1:3
    
    % Sum the images together
    tmp = nansum(dpn.d1.sp(:,:,mask(:,n) == 1),3);
    
    % Log / suppress
    if doLog
        tmp = log(tmp + 1);
    end
    if doSup && max(tmp(:)) > 0
        val = prctile(tmp(tmp > 0),95);
        if isnan(val)
            val = 0;
        end
        if val > 0 
            intmask = tmp > val;
            tmp(intmask) = val;
        end
    end
    
    % Update
    new(:,:,n) = real(tmp);
    
end

% Set the current axes and delete previous stuff
axes(fig.ax.mv);
f0 = get(fig.ax.mv,'Children');
delete(f0);

% Scale the image for plotting purposes
new = imScale(new);

% Here actually plot the things
imagesc(new);
set(gca,'YDir','reverse',...
    'XTick',[],...
    'YTick',[],...
    'XTickLabel',[],...
    'YTickLabel',[],...
    'LineWidth',5,...
    'TickLength',[0 0]);
xlim([0.5 sz(2)+0.5]);
ylim([0.5 sz(1)+0.5]);
title(fig.ax.mv,'RGB image');

% Change the other axes too
f0 = get(fig.ax.sp,'Children');
delete(f0);
title(fig.ax.sp,'');
set(fig.ax.sp,'XTick',[]);

% But only proceed if we have some annotations
if ~isfield(dpn,'anno')
    return
end

% Can we extract the annotated regions from this image?
[mask2,histID,~,] = desiAnnotationExtract(dpn);

% Reshape all variables of interest
tmp = reshape(new,[sz(1)*sz(2) 3]);
tmp = tmp(mask2 > 0,:);
histID = histID(mask2 > 0,:);

% Determine the colours
[~,colOrder] = sort(dpn.anno(:,5));
cols = vertcat(dpn.anno{:,3});
cols = cols(colOrder,:);

% Quick and dirty fix for multiple regions of the same annotation
[~,b,~] = unique(cols,'rows');
b = sort(b);
cols = cols(b,:);

if max(tmp(:)) > 0
    jsmBoxPlotMulti(tmp,histID,[1 2 3],...
        'Orientation',fig.ax.sp,...
        'Colours',cols,...
        'Legend',false,...%'Labels',{'Background','Blood','Muscle','Tumour'},...
        'Order',1:size(cols,1));

    set(fig.ax.sp,...
        'YTickLabel',{'R ','G ','B '},...
        'FontSize',10,...
        'FontWeight','normal',...
        'FontName','Helvetica');
end
xlabel(fig.ax.sp,'');



end