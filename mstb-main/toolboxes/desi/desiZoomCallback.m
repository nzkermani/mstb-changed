function desiZoomCallback(src,~,fig,value)
%desiZoomCallback - enable/disable the zoom in/out function

% Get the guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    set(src,'State','off');
    return
end

% Sort out the click boxes and return the handle to the clicked one
[stc] = clickBoxSort(fig,value);

% Now we need to decide which of the three outcomes is desired...
if strcmp(stc,'on') && value > 0
    
    % Zoom in
    disp('in');
    
    % Which axes to zoom into?
    hZoom = zoom(dpn.fig.fig);%,'Direction','in');
    
    % Set post-zoom callback
    set(hZoom,'ActionPostCallback',{@desiZoomPerform,dpn},...
        'ActionPreCallback',@desiZoomPre);

    % Enable the zoom
    set(hZoom,'Enable','on','Direction','in')
    
elseif strcmp(stc,'on') && value < 0
    
    % Zoom out
    disp('out');
    
    % Which axes?
    hZoom = zoom(dpn.fig.fig);%,'Direction','out');
    
    % Set post-zoom callback
    set(hZoom,'ActionPostCallback',{@desiZoomPerform,dpn},...
        'ActionPreCallback',@desiZoomPre);

    % Enable the zoom
    set(hZoom,'Enable','on','Direction','out')
    
else
    
    % Do nothing
    disp('stick');
    
    %zoom off;
    %hZoom = zoom(dpn.fig.fig);
    
    % Pre and post action callbacks
    %set(hZoom,'ActionPostCallback',[],'ActionPreCallback',[]);
    
    % Define the hZoom parameter, either that saved from the previous
    % in/out operation or the 
    try
        hZoom = dpn.fig.hZoom;
    catch
        hZoom = zoom(dpn.fig.fig);
    end
    
    
    % THis needs to be tested in the long term to ensure that there are no
    % problems with using it.
    try
        set(hZoom,'Enable','off');
    catch %err
        %err
        try
            zoom off
        catch %err2
            %err2
            zoom
        end
    end
    
end
    
% Add the hZoom to the dpn structure...
dpn.fig.hZoom = hZoom;
guidata(dpn.fig.fig,dpn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stc] = clickBoxSort(fig,value)
% Sort out the click boxes

% Which one was clicked?
if value > 0
    click = fig.tb.zoom(1);
    other = fig.tb.zoom(2);
else
    click = fig.tb.zoom(2);
    other = fig.tb.zoom(1);
end

% Get state of clicked box
stc = get(click,'State');

% If clicked button's state is 'on', then turn other to 'off'
if strcmp(stc,'on')
    set(other,'State','off');
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

