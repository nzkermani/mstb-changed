function fig = copyCurrentAxes(view)
%copyCurrentAxes


% Copy into new figure
copyobj(gca,figure);

% Get the axes handle
g = findobj('Parent',gcf,'Type','axes');

if nargin == 0
    view = 'full';
end

% Set the new position
switch view
    case 'full'
        set(g,'Position',[0 0 1 1]);
        
    case 'normal'
        set(g,'Position',[0.15 0.15 0.7 0.7]);
        
end



end

