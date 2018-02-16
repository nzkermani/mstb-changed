function pdSegmentDetermine(src,event,fig)
%pdSegmentDetermine - determine and show the segments on the image

% Get the data
data = guidata(fig.fig);
if isempty(data)
    return
end

% Now how many segments?
numS = str2double(fig.numSect.String);
cols = parula(numS);

% Let's delete any existing patches and outlines before we continue...
if isfield(data,'table')
    table = data.table;
    delete([table{:,1}]);
    delete([table{:,2}]);
    table = cell(numS,7);
else
    table = cell(numS,7);
end

% Anonymous function for coloured cell
colorgen = @(color,text) ['<html><table border=0 width=400 bgcolor=',...
    color,'><TR><TD>',text,'</TD></TR> </table></html>'];


% Get the image
img = fig.ax(2).CData;
img = img - min(img(:));
img = img / max(img(:));


% Now find them...
[ols] = pdSegmentOutlines(img,numS);
numS = size(ols,1);

% Resize the table in case there aren't enough outlines
table = table(1:numS,:);

% Now we need to show them on the image... And also to put them in the
% table and perhaps add rectangles as well...

% Let's first of all just draw them...
axes(fig.ax(1));
hold on;
for n = 1:numS
    
    % THis outline...
    ol = ols{n};
    
    % Determine the limits of the outlined areas...
    xl = [min(ol(:,2)) max(ol(:,2))] + [-5 5];
    yl = [min(ol(:,1)) max(ol(:,1))] + [-5 5];
    
    % Check within image boundaries
    if xl(1) < 1
        xl(1) = 1;
    end
    if yl(1) < 1;
        yl(1) = 1;
    end
    if xl(2) > size(img,2)
        xl(2) = size(img,2);
    end
    if yl(2) > size(img,1)
        yl(2) = size(img,1);
    end    
    
    % Now add a rectangle around the thing
    table{n,2} = patch([xl(1) xl(2) xl(2) xl(1)],...
        [yl(1) yl(1) yl(2) yl(2)],'green',...
        'FaceColor',cols(n,:),...
        'FaceAlpha',0.2,...
        'EdgeColor',cols(n,:));
    
    % Draw the outlines
    table{n,1} = plot(ol(:,2),ol(:,1),...
        'Color',cols(n,:),...
        'LineWidth',2);

    % What to put in here?
    table{n,4} = colorgen(rgb2hex(cols(n,:)),int2str(n));

    
    % Add information to the table
    table{n,3} = false;
    %table{n,4} = [1 1 1];
    table{n,5} = ['Section ' int2str(n)];
    
    % Save the boundaries
    table{n,6} = xl;
    table{n,7} = yl;
    
end

% Save the table...
data.table = table;

% Update the table...
fig.segTable.Data = table(:,3:5);

% Save the handles for the outlines and patches...
guidata(fig.fig,data);
    
    


end

