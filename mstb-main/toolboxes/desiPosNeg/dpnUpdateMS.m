function dpnUpdateMS(~,~,fig,flag)
%dpnUpdateMS - change the MS images in the window

% Get the guidata
dpn = guidata(fig.fig);

% Quit if there is nothing...
if isempty(dpn)
    return
end

% Ask the user if they are sure?
if nargin == 3
    decision = questdlg('Reset to beginning.  All annotations will be lost',...
        'Reset','Yes','No','No');
    if strcmp(decision,'No')
        return
    end
elseif nargin == 4
    % Then we force the refresh
end

% May want to remove the annotations by force as well...
if isfield(dpn,'anno')
    
    f0 = findobj('Type','patch');
    delete(f0);
    
    dpn = rmfield(dpn,'anno');
end

% Find enabled side menu buttons and disable them...
f0 = findobj('Tag','desiSideMenu');
set(f0,'State','off')
f0 = findobj('Parent',dpn.fig.pan2);
delete(f0);

% Restore and update d1
dpn.d1.img = log2(nansum(dpn.d1.sp,3)+1);
dpnIonImage([],[],fig.ax.ms1,dpn.d1.img);

% Do same for second image if we have it
if isfield(dpn,'d2')
    dpn.d2.img = nansum(dpn.d2.sp,3);
    dpnIonImage([],[],fig.ax.ms2,dpn.d2.img);
end

% Save the guidata...
guidata(fig.fig,dpn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

