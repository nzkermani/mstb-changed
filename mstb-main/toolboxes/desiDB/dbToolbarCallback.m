function dbToolbarCallback(~,~,fig,defP)
%dbCallbacks - where all of the callbacks for the database window are
%applied. Now a separate function so that the redraw function can avoid
%having to make a new window

% This is the button for the selection of all or none of the files
set(fig.selAll,'Callback', {@dbUIselAll,fig,[]});

% Change the path
set(fig.tb.new,'ClickedCallback', {@dbChangePath,fig});

% Set the proc button to draw the stuff in the panel
set(fig.tb.msaMake,...
    'ClickedCallback',{@dbPanelChange,fig,{@dbPanelMSA,fig,defP}},...
    'Tag','dbSideMenu');

% Now the button for inspecting a previous MSA
set(fig.tb.msaInspect,...
    'ClickedCallback',{@dbPanelChange,fig,{@dbPanelInspect,fig,defP}},...
    'Tag','dbSideMenu');

% Close all other windows
set(fig.tb.duck,'ClickedCallback',@closeOtherWindows);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
