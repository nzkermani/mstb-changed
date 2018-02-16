function dpnManipulateCallback(~,~,fig,mode,man)
%dpnManipulateCallback - change the ion image in the axes depending on the
%options specified in the mini window...

% Get the main figure's guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Determine ion mode
switch mode
    case 'pos'
        dx = dpn.d1;
        ax = fig.ax.ms1;
    case 'neg'
        dx = dpn.d2;
        ax = fig.ax.ms2;
end

% 1. ION IMAGE REFRESH...
ions = get(man.ion,'String');
ions = str2num(ions)'; %#ok<ST2NM>
ions = ions(:);
[dx.img,dx.opts.ions] = prepIonImage(dx,ions);

% In dual mode do we want to show the interpolated pixels?
if strcmp(dpn.mode,'dual')
    interpFlag = get(man.showInterp,'Value');
    
    if interpFlag
        dx.img(dx.isInterp) = 0;
    end
end

% Decide to log the image perhaps?
doLog = get(man.doLog,'Value');
if doLog
    os = nanmedian(dx.img(dx.img > 0));    
    dx.img = log(dx.img + os);    
end
dx.opts.doLog = doLog;

% Let's supress the high intensities
doSuppress = get(man.suppress,'Value');
if doSuppress
    maskVal = prctile(dx.img(:),95);
    mask = dx.img > maskVal;
    dx.img(mask) = maskVal;
end

% Finally update the ion image...
dpnIonImage([],[],ax,dx.img);

% Save dx back to the structure and then guidata, and add the annotations
switch mode
    case 'pos'
        dpn.d1 = dx;
        %dpnRedrawPatch([],[],fig,'ms1');
    case 'neg'
        dpn.d2 = dx;
        %dpnRedrawPatch([],[],fig,'ms2');
end

% Update fused ion image
if strcmp(dpn.mode,'dual')
    dpnUpdateFusedIonImage([],[],fig);
end

% Save the guidata
guidata(fig.fig,dpn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ii,ions] = prepIonImage(dx,ions)
% Prepare the ion image according to the ions selected...

numI = numel(ions);
if mod(numI,2) == 1
    ions = ions(1:end-1);
    numI = numI - 1;
end

mask = false(size(dx.mz));

for n = 1:2:numI
    
    tmp = dx.mz > ions(n) & dx.mz < ions(n+1);
    
    mask = mask + tmp;
    
end
mask = mask ~= 0;

ii = nansum(dx.sp(:,:,mask),3);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%