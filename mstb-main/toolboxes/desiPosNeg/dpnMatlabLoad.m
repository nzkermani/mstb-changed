function dpnMatlabLoad(imz,fig)
%dpnMatlabLoad - function that loads a saved file, and restores things as
%they should be...

% Open the file...
if isfield(imz,'dir')
    tmp = open([imz.dir imz.nam]);
    
    if isfield(tmp,'dpn')
        try        
            dpn = tmp.dpn;
        catch
            dpn = tmp;
        end
    elseif isfield(tmp,'data')
        dpn = tmp.data;
    end
    clear tmp;

else
    dpn = imz;
end

% If not, then we continue to load into the old fashioned window

% Ensure that older versions of .mat files are compatible where possibe...
if ~isfield(dpn,'mode')
    if isfield(dpn,'d2')
        dpn.mode = 'dual';
    else
        dpn.mode = 'single';
    end
end

% Update the file information to represent the imz file
if isfield(dpn,'file')
    imz.raw = dpn.file;
elseif isfield(dpn,'imz')
    imz.raw = dpn.imz;
end
%dpn.file = imz;

% Add to the guidata
dpn.fig = fig;
guidata(fig.fig,dpn);

% Change the layout
if strcmp(dpn.mode,'single')
    if strcmp(get(fig.tb.layout,'State'),'on')
        set(fig.tb.layout,'State','off');
        desiChangeLayout(fig.tb.layout,[],fig);
    end
end

% Set all menus to be off...
f1 = findobj('Tag','desiSideMenu');
set(f1,'State','off');

% Delete anything in the side menu panel
f0 = get(fig.pan2,'Children');
delete(f0);

% Load the optical image
coregFlag = true;
if ~isfield(dpn,'opt')
    try
        dirIcon = deSlash([pwd filesep 'desi/icons/']);
        load([dirIcon filesep 'Images.mat']);
        idx = ceil(rand(1,1) * numel(images)); %#ok<USENS>
        opt = images{1,idx};
    catch
        opt = rand(10,10);
    end
    coregFlag = false;    
else
    if isfield(dpn.opt,'coreg')
        opt = dpn.opt.coreg;    
    elseif isfield(dpn.opt,'orig')
        opt = dpn.opt.orig;
    else
        error('No optical image found');
    end
end

% Load the positive MS data
if isfield(dpn.d1,'img')
    ms1 = dpn.d1.img;
else
    try        
        ms1 = nansum(dpn.d1.sp,3);
        dpn.d1.img = ms1;
        guidata(fig.fig,dpn);
    catch
        error('No MS1');
    end
end

% Load the negative MS data
if strcmp(dpn.mode,'dual')
    if isfield(dpn.d2,'img')
        ms2 = dpn.d2.img;
    else
        error('No MS2');
    end
else
    ms2 = [];
end

% Draw the optical image...
set(fig.ax.opt(2),'CData',opt);
set(fig.ax.opt(1),...
    'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(opt,2)],...
    'YLim',[1 size(opt,1)],...
    'TickLength',[0 0]);

set(fig.ax.ms1(2),'CData',ms1);
set(fig.ax.ms1(1),...
    'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(ms1,2)],...
    'YLim',[1 size(ms1,1)],...
    'TickLength',[0 0]);

if strcmp(dpn.mode,'dual')
    set(fig.ax.ms2(2),'CData',ms2);
    set(fig.ax.ms2(1),...
        'XTick',[],...
        'YTick',[],...
        'XLim',[1 size(ms2,2)],...
        'YLim',[1 size(ms2,1)],...
        'TickLength',[0 0]);
end

% Update the fused ion image to end...
if strcmp(dpn.mode,'dual')
    dpnUpdateFusedIonImage([],[],fig);
end

colormap(redbluecmap);

if coregFlag
    set(dpn.fig.tb.coreg,'Enable','on');
    if strcmp(dpn.mode,'single')
        set(dpn.fig.tb.fiducial,'Enable','on');
    end
else
    set(dpn.fig.tb.coreg,'Enable','off');
    if strcmp(dpn.mode,'single')
        set(dpn.fig.tb.fiducial,'Enable','off');
    end
end

% Later on would come things like the annotations and stuff, but I need to
% do something about that first!
dpnRedrawPatch([],[],fig,'opt');
dpnRedrawPatch([],[],fig,'ms1');
if strcmp(dpn.mode,'dual')
    dpnRedrawPatch([],[],fig,'ms2');
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = redrawPatch(x,y,col)
% Draw a patch in the CURRENT axes

if numel(x) == 1 && numel(y) == 1
    % Then this is a single scatter point of annotation
    h = scatter(x(1),y(1),200,col,'o','filled');
    
else
    % Then multiple points, so we patch them together
    h = patch(x,y,...
        col,...
        'EdgeColor',col,...
        'FaceColor',col,...
        'FaceAlpha',0.4,...
        'LineWidth',3);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

