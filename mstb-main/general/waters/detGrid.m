%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [crx,cry] = detGrid(xy)
% From the locations of points, create a regular grid

minX = min(xy(:,1));
minY = min(xy(:,2));

maxX = max(xy(:,1));
maxY = max(xy(:,2));

resX = mode(diff(unique(xy(:,1))));
resY = mode(diff(unique(xy(:,2))));

% Ensure that they are rounded to just the one decimal place
resX = str2double(sprintf('%0.1f',resX));
resY = str2double(sprintf('%0.1f',resY));

% Convert the xy locations to match those in gx and gy
crx = round((xy(:,1)-minX) /resX) + 1;
cry = round((xy(:,2)-minY) /resY) + 1;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
