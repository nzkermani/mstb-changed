function desiCrop(~,~,fig)
%desiCrop - get the user to specify the region that they want to keep.  All
%other parts will be removed permanently.

% Define which method to use: either rect (@getrect) or gin (@ginput)
method = 'rect'; 

% Guidata gather
dpn = guidata(fig.fig);
if isempty(dpn)
    return;
end

% Ensure also that we have imported some MS data
if ~isfield(dpn,'d1')
    return;
end

% Add a large-ish warning dialog box to inform the user that this will a
% permanent change to the file, which cannot be undone
q1 = ['CROP THE MS IMAGE' char(10) 'This function allows you to '...
    'crop the MS image by drawning a rectangle around the region '...
    'that you would like to keep.  It is NOT reversible without '...
    're-processing from the imzML/raw data.  Do you wish to proceed?'];

choice = questdlg(q1,...
    'Cropped Happiness',...
    'No','Yes','No');

switch choice
    case 'No'
        return
    otherwise
        % Continue
end


% Determine the actual sizeo f the image
sz = size(dpn.d1.img);

% Now that we have some data, set the axes to be the MS1 axes and then use
% the ginput function to get the region

% Set the axes
axes(dpn.fig.ax.ms1(1));

% Region definition
switch method
    
    case 'gin'
               
        [x,y] = ginput(2);
        if x(1) < 1
            x(1) = 1;
        end
        if y(1) < 1
            y(1) = 1;
        end
        
    case 'rect'

        [rr] = getrect;
        
        if rr(1) < 1
            rr(1) = 1;
        end
        if rr(2) < 1
            rr(2) = 1;
        end
                       
        % Convert these to xx and yy coordinates
        x = [rr(1) rr(3) + rr(1)];
        y = [rr(2) rr(4) + rr(2)];
        
end

% Just round to the nearest pixel
x = round(x);
y = round(y);

% Now ensure that the coordinates don't exceed the width of the image
if x(2) > sz(2)
    x(2) = sz(2);
end
if y(2) > sz(1)
    y(2) = sz(1);
end

% x and y now finally are useful. We should draw a patch on the axes an
% then ask the user to confirm that they are happy with the region shown.
% If yes, then they can click continue / or not etc

% Patch
px = [x(1) x(1) x(2) x(2)];
py = [y(1) y(2) y(2) y(1)];

col = [0.5 0.5 0.5];
p = patch(px,py,col,...
    'EdgeColor',col,...
    'FaceColor',col,...
    'FaceAlpha',0.4,...
    'LineWidth',3);

% Ask the user to confirm happiness
choice = questdlg('Are you happy?',...
    'Cropped Happiness',...
    'No','Yes','No');

switch choice
    
    case 'Yes'
        
        dpn.d1.sp = dpn.d1.sp(y(1):y(2),x(1):x(2),:);
        dpn.d1.img = dpn.d1.img(y(1):y(2),x(1):x(2));
        %dpn.d1.tobg = dpn.d1.tobg(y(1):y(2),x(1):x(2));
        
        
        guidata(dpn.fig.fig,dpn);
        
        % Delete the patch that we drew
        delete(p);
        
        % Call the reset function...
        dpnUpdateMS([],[],dpn.fig,true);
        
    case 'No'
        delete(p);
end


end

