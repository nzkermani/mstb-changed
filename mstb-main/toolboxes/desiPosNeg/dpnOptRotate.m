function dpnOptRotate(~,~,fig,cor,direction)
%dpnOptRotate - get the optical image for the giudata, rotate it according
%to the instruction, and redraw the image

% Get guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Get the optical image
opt = dpn.opt.orig;
opt = rot90(opt,direction);
%opt = imrotate(opt,direction);

% Save the optical image to the structure
dpn.opt.orig = opt;
guidata(fig.fig,dpn);

% Now that we have an image then plot it...
set(fig.ax.opt(2),'CData',dpn.opt.orig,...
    'XData',[1 size(dpn.opt.orig,2)],...
    'YData',[1 size(dpn.opt.orig,1)]);

% Change the axes properties
set(fig.ax.opt(1),...
    'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(dpn.opt.orig,2)],...
    'YLim',[1 size(dpn.opt.orig,1)]);

% I'm not sure that I know what this actually is for?
try
    set(cor.listOP,'Value',1);
catch
    disp('could not set');
end

end

