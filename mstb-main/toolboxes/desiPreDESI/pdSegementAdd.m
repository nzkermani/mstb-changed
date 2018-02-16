function pdSegementAdd(~,~,fig)
%pdSegmentAdd - allow the user to draw a new segment

% Get guidata
data = guidata(fig.fig);
if isempty(data)
    return
end

% Determine which colour it could be...
origSeg = str2double(fig.numSect.String);
cols = parula(origSeg);

% Labels...
nos = zeros(size(data.table,1),1);
for n = 1:numel(nos)
    
    f1 = strfind(data.table{n,4},'<TD>');
    f2 = strfind(data.table{n,4},'</TD>');
    
    txt = data.table{n,4}(f1+4:f2-1);
    
    nos(n,1) = str2double(txt);

end
nos = sort(nos);
% So what is the lowest integer value that remains?
if numel(nos) >= origSeg
    newSeg = numel(nos) + 1;
elseif numel(nos) < origSeg
    tmp = 1:max(nos);
    newSeg = setdiff(tmp,nos);
    if numel(newSeg) == 0
        newSeg = numel(nos) + 1;
    end
elseif size(nos,1) == 0
    newSeg = 1;
end

% This will be its number
newSeg = newSeg(1);

% What will its colour be?
if newSeg <= size(cols,1)
    thisCol = cols(newSeg,:);
else
    thisCol = rand(1,3);
end

% Table index...
ti = size(data.table,1) + 1;
    
% Now we need to ask the user to draw a rectangle
rect = getrect(fig.ax(1));

% Will have to round off the values for whole integer pixels
rect = round(rect);

% How about drawing it...
xx = [rect(1) rect(1)+rect(3) rect(1)+rect(3) rect(1)];
yy = [rect(2) rect(2) rect(2)+rect(4) rect(2)+rect(4)];

% Ensure that it doesn't stray
img = fig.ax(2).CData;
if xx(1) < 1
    xx([1 4]) = 1;
end
if yy(1) < 1
    yy([1 2]) = 1;
end
if xx(2) > size(img,2);
    xx([2 3]) = size(img,2);
end
if yy(3) > size(img,1)
    yy([3 4]) = size(img,1);
end

data.table{ti,2} = patch(xx,yy,'green',...
        'FaceColor',thisCol,...
        'FaceAlpha',0.2,...
        'EdgeColor',thisCol);

% Anonymous function for coloured cell
colorgen = @(color,text) ['<html><table border=0 width=400 bgcolor=',...
    color,'><TR><TD>',text,'</TD></TR> </table></html>'];

% Add the colour information
data.table{ti,4} = colorgen(rgb2hex(thisCol),int2str(newSeg));
   
% Add information to the table
data.table{ti,3} = false;
data.table{ti,5} = ['Section ' int2str(newSeg)];

% Save the boundaries
data.table{ti,6} = xx([1 2]);
data.table{ti,7} = yy([1 3]);

% Update the table in the gui
fig.segTable.Data = data.table(:,3:5);

% Update guidata
guidata(fig.fig,data);

end

