function dbPanelChange(src,~,fig,fcn)
%dbPanelChange - flick between variously triggered panels. Only gets
%complicated with various co-existing functions, such as in the desi
%toolbox. Currently, this is sufficiently simple.

% Get the state of the button
state = get(src,'State');
        
% Get current children of the panel
f0 = get(fig.panel,'Children');
delete(f0);

% Just return if the button has been turned off
if strcmp(state,'off')    
    return
end

% Turn all other (relevant) buttons off
f1 = findobj('Tag','dbSideMenu');
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

