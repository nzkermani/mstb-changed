function desiZoomReset(~,~,fig)
%desiZoomReset - restore both axes to show all of the images

% Guidata gather
dpn = guidata(fig.fig);
if isempty(dpn)
    return;
end

% Set both zoom buttons to be off
set(dpn.fig.tb.zoom,'State','off','Enable','on');
desiZoomCallback([],[],dpn.fig,0);

% Change the view on the optical image
if isfield(dpn.opt,'coreg')
    set(dpn.fig.ax.opt(1),...
        'XLim',[1 size(dpn.opt.coreg,2)],...
        'YLim',[1 size(dpn.opt.coreg,1)]);
    %grid off;
end

% This changes the axes
set(dpn.fig.ax.ms1(1),...
    'XLim',[1 size(dpn.d1.img,2)],...
    'YLim',[1 size(dpn.d1.img,1)]);
%grid off;

% Also potentially the other one
if strcmp(dpn.mode,'dual')
    set(dpn.fig.ax.ms2(1),...
        'XLim',[1 size(dpn.d2.img,2)],...
        'YLim',[1 size(dpn.d2.img,1)]);
    
    set(dpn.fig.ax.fu(1),...
        'XLim',[1 size(dpn.d2.img,2)],...
        'YLim',[1 size(dpn.d2.img,1)]);
    
end    


% Also zoom out of the other axes if we can...
if strcmp(dpn.mode,'single')
    if strcmp(get(fig.tb.layout,'State'),'on') 
        zoomOutSecondaryAxes(dpn.fig.ax.mv,dpn);
        zoomOutSecondaryAxes(dpn.fig.ax.sp,dpn);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zoomOutSecondaryAxes(parent,dpn)

f0 = get(parent,'Children');

% This is an image...
f1 = find(strcmp(get(f0,'Type'),'image'));
if numel(f1) > 0
    sz = get(f0(f1),'CData');
    sz = size(sz);

    % See if this is the same size as the actual image
    if sum(sz(1:2) == size(dpn.d1.img)) == 2

        % MS image sized
        xlim(parent,[1 sz(2)]);
        ylim(parent,[1 sz(1)]);

    else

        % Heat map sized
        xlim(parent,[0.5 sz(2)+0]);
        ylim(parent,[0 sz(1)+1]);
    end

else
    % Then there is no image here, so we need to consider that it is
    % loadings and stuff...
    xlim(parent,'auto');
    ylim(parent,'auto');
    
    %get(f0,'Type')
    %disp('zoom out please');
end
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

