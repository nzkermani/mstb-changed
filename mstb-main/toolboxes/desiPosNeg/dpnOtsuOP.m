function dpnOtsuOP(~,~,fig,cor)
%dpnOtsuOP - do Otsu thresholding on the OP image

% Get the guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Determine which mode is to be thresholded
tmp = get(cor.listOP,'String');
mode = get(cor.listOP,'Value');
mode = tmp{mode};
switch mode
    case 'Original'
        % Just need to ensure that we draw the optical image
        set(dpn.fig.ax.opt(2),'CData',dpn.opt.orig,...
            'XData',[1 size(dpn.opt.orig,2)],...
            'YData',[1 size(dpn.opt.orig,1)]);
        set(dpn.fig.ax.opt(1),'XTick',[],...
            'YTick',[],...
            'XLim',[1 size(dpn.opt.orig,2)],...
            'YLim',[1 size(dpn.opt.orig,1)]);
        return
        
    case 'Otsu'
        % Just continue...
end

% The optical image will need resizing from the original, down to the size
% of the MS image
lowRes = zeros([size(dpn.d1.img) 3]);
for n = 1:3
    lowRes(:,:,n) = imresize(dpn.opt.orig(:,:,n),size(dpn.d1.img));
end
dpn.opt.lowRes = lowRes ./ max(lowRes(:));

% Convert this to grayscale
dpn.opt.gray = nansum(lowRes,3);
dpn.opt.gray = max(dpn.opt.gray(:)) - dpn.opt.gray;

% What is the state of the slider at the moment?
val = round(get(cor.slideOP,'Value'));

% Run the Otsu tresholding method
[dpn.opt.tobg,~] = dpnTOBG(dpn.opt.gray,[],val);

% Perhaps we can consider smoothing the image here...
filt = fspecial('average',3);
dpn.opt.tobg = filter2(filt,dpn.opt.tobg) >  0.75;

% Save to the guitata
guidata(fig.fig,dpn);

% Draw the things in the axes...
set(dpn.fig.ax.opt(2),'CData',dpn.opt.tobg,...
    'XData',[1 size(dpn.opt.tobg,2)],...
    'YData',[1 size(dpn.opt.tobg,2)]);
set(dpn.fig.ax.opt(1),'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(dpn.opt.tobg,2)],...
    'YLim',[1 size(dpn.opt.tobg,1)]);

end

