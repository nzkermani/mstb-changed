function desiDrawSideMenu(src,~,fig,fcn)
%desiSideMenu - draw or clear the side menu, and run the appropriate
%function.

% Exception for loading a new file
if strcmp(char(fcn{1}),'desiFileUpload') || strcmp(char(fcn{1}),'desiQQQUpload')
    newfile = true;
else
    newfile = false;
end

% Guidata
dpn = guidata(fig.fig);
if isempty(dpn) && ~newfile
    set(src,'State','off');
    return
end

% Want to restore the ion images when leaving the side menus.  Other axes
% also could do with reseting too
if ~newfile
    dpnIonImage([],[],fig.ax.ms1,dpn.d1.img);
    if strcmp(dpn.mode,'dual')
        dpnIonImage([],[],fig.ax.ms2,dpn.d2.img);
        dpnUpdateFusedIonImage([],[],fig);
    end
end

% What about deleting any fiducial markers that may exist...
f0 = findobj('Tag','fiducialMarker');
delete(f0);

% Reset the axes
desiAxesReset(fig.ax.mv);
desiAxesReset(fig.ax.sp);

% Get the state of the button
state = get(src,'State');

% If off, then remove items from the side menu (fig.pan2).
% If on, then we need to turn off all other side-menu buttons...
% ...and then draw/run the desired function
        
% Actually we need to do this for all functions whether they are on or off.
% The delete part has been moved after this 'if'
f0 = get(fig.pan2,'Children');

% Return if 'off' once we have delete the residue
if strcmp(state,'off') && strcmp(char(fcn{1}),'dpnAnnotate')
    
    % Here we want to run the function that updates the annotation table,
    % as now we don't have the close request function or the 'Finish'
    % button
    man = get(src,'UserData');
    dpnFinalAnnotation(src,[],fig,man);
    
    % The delete must come after the above function
    delete(f0);
    return
    
elseif strcmp(state,'off') && strcmp(char(fcn{1}),'dpnCoreg')
    
    % Here we must ensure that the images are changed to show the ion and
    % optical images, rather than being stuck with the Otsu thresholded
    % images
    try
        set(fig.ax.opt(2),'CData',dpn.opt.coreg);    
        set(fig.ax.opt(1),'XTick',[],...
            'YTick',[],...
            'XLim',[1 size(dpn.opt.coreg,2)],...
            'YLim',[1 size(dpn.opt.coreg,1)]);

    catch
        set(fig.ax.opt(2),'CData',dpn.opt.orig);
        set(fig.ax.opt(1),'XTick',[],...
            'YTick',[],...
            'XLim',[1 size(dpn.opt.orig,2)],...
            'YLim',[1 size(dpn.opt.orig,1)]);

        disp('No coreg performed');
    end
    dpnIonImage([],[],fig.ax.ms1,dpn.d1.img);

    
    delete(f0);
    return
    
elseif strcmp(state,'off') && strcmp(char(fcn{1}),'dpnSegmentationWindow')
    
    % Here we need to return to the original ion images if we are in dual
    % mode. Note sure about single ion mode as it does use the same
    % function
    if strcmp(dpn.mode,'dual')
        
        dpnIonImage([],[],fig.ax.ms1,dpn.d1.img);
        dpnIonImage([],[],fig.ax.ms2,dpn.d2.img);
        dpnUpdateFusedIonImage([],[],fig);
    else
        
    end
    
    delete(f0);
    return
   
elseif strcmp(state,'off')
    
    % Replicated code as delete can't happen too soon, and the function
    % must return from here
    delete(f0);
    return
else
    
    % Delete, but no return
    delete(f0);
end

% Turn all other buttons off
f1 = findobj('Tag','desiSideMenu');
f1 = setdiff(f1,src);
set(f1,'State','off');

% If 'on', then we need to just run the function as specified by the input
% 'fcn'.  There must be a better way to run this!
numA = numel(fcn);
switch numA
    case 1
        feval(fcn{1},src,[]);
    case 2
        feval(fcn{1},src,[],fcn{2});
    case 3
        feval(fcn{1},src,[],fcn{2},fcn{3});
    case 4
        feval(fcn{1},src,[],fcn{2},fcn{3},fcn{4});
    case 5
        feval(fcn{1},src,[],fcn{2},fcn{3},fcn{4},fcn{5});
    otherwise
        error('The function is too long');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

