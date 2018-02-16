function xxxRotateMSI(~,~,fig,cor,direction)
%xxxRotateMSI - get the optical image for the giudata, rotate it according
%to the instruction, and redraw the image

% Get guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Change axes view
xxxEnforceLayoutChange(fig,'single','off');

% Clear the axes... but what if we have annotation boxes?
f0 = findobj('Type','patch');
delete(f0);
f0 = findobj('Type','scatter');
delete(f0);
dpn.anno = [];

% Remove existing grid
try
    dpn.fig = rmfield(dpn.fig,'grid');
catch
end

% Reset in case otsu threshold has been selected
cor.listMS.Value = 1;

% Rotate dpn.d1
dpn.d1.sp = rot90(dpn.d1.sp,direction);

if isfield(dpn.d1,'img')
    dpn.d1.img = rot90(dpn.d1.img,direction);
end
if isfield(dpn.d1,'tobg')
    dpn.d1.tobg = rot90(dpn.d1.tobg,direction);
end

% Rotate optical image formats
if isfield(dpn.opt,'lowRes')
    dpn.opt.lowRes = rot90(dpn.opt.lowRes,direction);
end
if isfield(dpn.opt,'gray')
    dpn.opt.gray = rot90(dpn.opt.gray,direction);
end
if isfield(dpn.opt,'tobg')
    dpn.opt.tobg = rot90(dpn.opt.tobg,direction);
end
if isfield(dpn.opt,'coreg')
    dpn.opt.coreg = rot90(dpn.opt.coreg,direction);
end

% Save it
guidata(fig.fig,dpn);

% Update the ion image
dpnIonImage([],[],fig.ax.ms1,dpn.d1.img);

% Now that we have an image then plot it...
if isfield(dpn.opt,'coreg')
    set(fig.ax.opt(2),'CData',dpn.opt.coreg,...
        'XData',[1 size(dpn.opt.coreg,2)],...
        'YData',[1 size(dpn.opt.coreg,1)]);
    
    % Change the axes properties
    set(fig.ax.opt(1),...
        'XTick',[],...
        'YTick',[],...
        'XLim',[1 size(dpn.opt.coreg,2)],...
        'YLim',[1 size(dpn.opt.coreg,1)]);
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

