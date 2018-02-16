function desiAxesReset(ax)
%desiAxesReset - find and delete children in the supplmentary 'detail' axes

f0 = get(ax,'Children');
delete(f0);
title(ax,'');
set(ax,'XTick',[],...
    'YTick',[],...
    'LineWidth',5);

end

