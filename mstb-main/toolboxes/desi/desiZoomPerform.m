function desiZoomPerform(~,axes,dpn)
%desiZoomPerform - change the figure according to the zoom operation...

%disp('DESI Zoom Callback');

% Get the axes limits
xl = get(axes.Axes,'XLim');
yl = get(axes.Axes,'YLim');

% Is this dual or single ion mode?
if strcmp(dpn.mode,'dual')
    flag = true;
    [ax,sz,zoomFlag] = zoomDualMode(axes,dpn);
else
    flag = false;
    [ax,sz,zoomFlag] = zoomSingleMode(axes,dpn);
end

% Simple check in case it all fails badly
if ~zoomFlag
    disp('No special zoom mode');
    return
    %error('Zoom issue here!');
end


% Convert to % of main axes
xlp = xl / sz(1,2);
ylp = yl / sz(1,1);

% Convert to pixels of the other axes
newx = xlp .* sz(2,2);
newy = ylp .* sz(2,1);

% Ensure that we haven't gone over the limits of the axes
if newx(1) < 0.5
    newx(1) = 0.5;
end
if newy(1) < 0.5
    newy(1) = 0.5;
end
if newx(2) == sz(2,2)
    newx(2) = newx(2) + 0.5;
end
if newy(2) == sz(2,1)
    newy(2) = newy(2) + 0.5;
end


switch ax{2}
    case 'MS'
        xlim(dpn.fig.ax.ms1(1),newx);
        ylim(dpn.fig.ax.ms1(1),newy);
        
        % Do this for the second ion image
        if flag
            
            % MS2
            xlim(dpn.fig.ax.ms2(1),newx);
            ylim(dpn.fig.ax.ms2(1),newy);
            
            % FU
            xlim(dpn.fig.ax.fu(1),newx);
            ylim(dpn.fig.ax.fu(1),newy);
        end
            
        
    case 'Opt'
        xlim(dpn.fig.ax.opt(1),newx);
        ylim(dpn.fig.ax.opt(1),newy);
        
end

% This used to go in the space after the zoomFlag check
% % Determine which axes was zoomed in/out of
% if axes.Axes == dpn.fig.ax.sp
%     ax = {'Spectra'};
%     sz = [];
%     return
%     
% elseif axes.Axes == dpn.fig.ax.mv
%     % Maybe could determine if this is an image too...
%     
%     % If it was an image we look a little harder    
%     f0 = get(axes.Axes,'Children');
%     if numel(f0) == 1 && strcmp(get(f0,'Type'),'image')
%         ax = {'MS','Opt'};
%         sz = [size(dpn.d1.sp); size(dpn.opt.coreg)];
%         ex = true;
%         
%         disp('is an image');
%         %return
%     
%     else
%         ax = {'MV'};
%         sz = [];
%         ex = false;
%         return
%     end
%     return
%     
% elseif axes.Axes == dpn.fig.ax.opt(1)
%     ax = {'Opt','MS'};
%     sz = [size(dpn.opt.coreg); size(dpn.d1.sp)];
% elseif axes.Axes == dpn.fig.ax.ms1(1)
%     ax = {'MS','Opt'};
%     sz = [size(dpn.d1.sp); size(dpn.opt.coreg)];
% elseif axes.Axes == dpn.fig.ax.ms2(1) % fail in single mode here
%     ax = {'MS','Opt'};
%     sz = [size(dpn.d2.sp); size(dpn.opt.coreg)];
% else
%     disp('some other axes');
%     return
% end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ax,sz,zoomFlag] = zoomSingleMode(axes,dpn)

zoomFlag = true;

% Determine which axes was zoomed in/out of
if axes.Axes == dpn.fig.ax.sp
    ax = {'Spectra'};
    sz = [];
    zoomFlag = false;
    %return
    
elseif axes.Axes == dpn.fig.ax.mv
    % Maybe could determine if this is an image too...
    
    % If it was an image we look a little harder    
    f0 = get(axes.Axes,'Children');
    if numel(f0) == 1 && strcmp(get(f0,'Type'),'image')
        ax = {'MS','Opt'};
        sz = [size(dpn.d1.sp); size(dpn.opt.coreg)];
        ex = true;
        
        disp('is an image');
        %return
    
    else
        ax = {'MV'};
        sz = [];
        ex = false;
        %zoomFlag = false;
        %return
    end
    zoomFlag = false;
    %return
    
elseif axes.Axes == dpn.fig.ax.opt(1)
    ax = {'Opt','MS'};
    sz = [size(dpn.opt.coreg); size(dpn.d1.sp)];
elseif axes.Axes == dpn.fig.ax.ms1(1)
    ax = {'MS','Opt'};
    sz = [size(dpn.d1.sp); size(dpn.opt.coreg)];
elseif axes.Axes == dpn.fig.ax.ms2(1) % fail in single mode here
    ax = {'MS','Opt'};
    sz = [size(dpn.d2.sp); size(dpn.opt.coreg)];
else
    disp('some other axes');
    ax = {'Unknown'};
    sz = [];
    zoomFlag = false;    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ax,sz,zoomFlag] = zoomDualMode(axes,dpn)

zoomFlag = true;

% Determine which axes was zoomed in/out of
if axes.Axes == dpn.fig.ax.opt(1)
    ax = {'Opt','MS'};
    sz = [size(dpn.opt.coreg); size(dpn.d1.sp)];
elseif axes.Axes == dpn.fig.ax.ms1(1)
    ax = {'MS','Opt'};
    sz = [size(dpn.d1.sp); size(dpn.opt.coreg)];
elseif axes.Axes == dpn.fig.ax.ms2(1) % fail in single mode here
    ax = {'MS','Opt'};
    sz = [size(dpn.d2.sp); size(dpn.opt.coreg)];
elseif axes.Axes == dpn.fig.ax.fu(1)
    ax = {'MS','Opt'};
    sz = [size(dpn.d1.sp); size(dpn.opt.coreg)];
else
    disp('some other axes');
    zoomFlag = false;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeAxes()


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
