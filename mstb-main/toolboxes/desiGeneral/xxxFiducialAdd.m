function xxxFiducialAdd(~,~,fig,man)
% xxxFiducualAdd - add a single point on each axes

% Get the guidata
dpn = guidata(fig.fig);

% Define the two parents in order...
par1 = dpn.fig.ax.opt(1);
par2 = dpn.fig.ax.ms1(1);

% Determine how many annotations have previously been made...
if isfield(dpn.opt,'fiducial')
    if isempty(dpn.opt.fiducial)
        isFirst = true;
        idx = 1;
    else
        isFirst = false;

        % Max annotations?    
        maxAnno = dpn.opt.fiducial{end,3};

        % Current annotations
        curAnno = cell2mat(dpn.opt.fiducial(:,3));

        % Differences between the two
        sdf = setdiff(1:maxAnno,curAnno);
        if isempty(sdf)
            idx = maxAnno + 1;
        else
            idx = sdf(1);
        end
    end    
else
    isFirst = true;
    idx = 1;
end

% Ensure that we can't have too many markers - there are only 15 colours
% and this many shouldn't be necessary...
if idx > 15
    disp('No more markers are allowed - delete others');
    return
end

% Ask the user to select points on each axes
axes(par1);
[x1,y1] = ginput(1);
hold on;
h1 = scatter(x1,y1,200,'o','MarkerFaceColor',man.cols(idx,:),...
    'MarkerEdgeColor','k',...
    'Tag','fiducialMarker');

axes(par2);
[x2,y2] = ginput(1);
hold on;
h2 = scatter(x2,y2,200,'d','MarkerFaceColor',man.cols(idx,:),...
    'MarkerEdgeColor','k',...
    'Tag','fiducialMarker');
       
% Put the coordinates in a cell
tabl = {h1 h2 x1 y1 x2 y2};

% Anonymous function for coloured cell
colorgen = @(color,text) ['<html><table border=0 width=400 bgcolor=',...
    color,'><TR><TD>',text,'</TD></TR> </table></html>'];
col2 = colorgen(rgb2hex(man.cols(idx,:)),int2str(idx));

% Now we need to add this to the table
line = [false col2 idx man.cols(idx,:) tabl];
if isFirst
    dpn.opt.fiducial = line;
else
    dpn.opt.fiducial = cat(1,dpn.opt.fiducial,line);
end

% Sort the table to account for deletions and various other things...
dpn.opt.fiducial = sortrows(dpn.opt.fiducial,3);

% Now change the table data
set(man.tab,'Data',dpn.opt.fiducial(:,1:2));

% Now consider updating the guidata
guidata(fig.fig,dpn);


end

