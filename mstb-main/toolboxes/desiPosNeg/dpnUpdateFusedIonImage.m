function dpnUpdateFusedIonImage(~,~,fig)
%dpnUpdateFusedImage - change the fused image to represent the ion images
%shown in the individual axes

% Get the main figure's guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Get individual axes
im1 = get(fig.ax.ms1(2),'CData');
im2 = get(fig.ax.ms2(2),'CData');

% Sum together
imF = im1/max(im1(:)) + im2/max(im2(:));

% Add to the axes
dpnIonImage([],[],fig.ax.fu,imF);

end

