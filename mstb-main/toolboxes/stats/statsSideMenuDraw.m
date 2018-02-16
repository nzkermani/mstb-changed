function statsSideMenuDraw(src,~,fig,fcn)
%statsSideMenuDraw - draw or clear the side menu, and run the appropriate
%function.

% Guidata
sts = guidata(fig.fig);
if isempty(sts)
    set(src,'State','off');
    return
end

% Get the state of the button
state = get(src,'State');

% If off, then remove items from the side menu (fig.pan2).
% If on, then we need to turn off all other side-menu buttons...
% ...and then draw/run the desired function
        
% Actually we need to do this for all functions whether they are on or off.
% The delete part has been moved after this 'if'
f0 = get(fig.pan2,'Children');

% Return if 'off' once we have deleted the residue
if strcmp(state,'off') && strcmp(char(fcn{1}),'statsLayoutChange')
    
    % Delete the side menu
    delete(f0);
    
    % Add in the correct layout function
    fcn{end+1} = sts.datatype;
    
    % Still need to run the function to visualise the various bits
    %statsLayoutChange(src,~,fig)
    
elseif strcmp(state,'on') && strcmp(char(fcn{1}),'statsLayoutChange')
    
    % Add in the correct layout function to help with layout challenges
    fcn{end+1} = sts.datatype;
    
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
f1 = findobj('Tag','statsSideMenu');
f1 = setdiff(f1,src);
set(f1,'State','off');

% If the table is still showing then remove it if another function is
% called. Perhaps there is a better way to deal with this, but here is the
% quick and dirty fix
if ~strcmp(char(fcn{1}),'statsLayoutChange') && strcmp(fig.tb.table.State,'on')
    fig.tb.table.State = 'off';
    statsLayoutChange(fig.tb.table,[],fig,sts.datatype);
end

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

