function dpnDrawAnnotation(src,~,fig,man)
%dpnDrawAnnotation - draw a thingy on the optical image and go from
%there...

% Get the square's colour
col = get(src,'BackgroundColor');
usd = get(src,'UserData');

% Guidata
dpn = guidata(fig.fig);

%Set all pushbuttons in the panel to be disabled, so only one annotation
%can be performed
f0 = findobj('Parent',dpn.fig.pan2,'Style','pushbutton');
set(f0,'Enable','off');

% Get a rectangle (or single point) on the optical image
method = get(man.rect,'Value');
if method == 1
    method = 'rect';
else
    method = 'gin';
end

% Fixed size?
fixed = get(man.fixed,'Value') == 1;
if fixed
    method = 'fixed';
    fixVal = 6;
end
    
% Ensure that a grid has been calculated - if not, then do it now
if ~isfield(dpn.fig,'grid')
    desiDetermineGrid([],[],fig);
    dpn = guidata(fig.fig);
end

% Apply one of the methods
switch method
    case 'rect'
        rect = getrect(dpn.fig.ax.opt(1));

        % Convert to locations
        xv = round([rect(1) rect(1)+rect(3)]);
        yv = round([rect(2) rect(2)+rect(4)]);

    case 'gin'
        axes(dpn.fig.ax.opt(1));
        [x1,y1] = ginput(2);
        
        % Convert the values in rect to locations
        xv = sort(round(x1'));
        yv = sort(round(y1'));
        
    case 'fixed'
        
        % Note that this method needs to differentiate as to whether the
        % optical image is MS (same size) or H&E (diff size)
        if size(dpn.d1.img,1) == size(dpn.opt.coreg,1)
            optMode = 'sameSize';
        else
            optMode = 'diffSize';
        end
        
        % Set the axes as the optical image
        axes(dpn.fig.ax.opt(1));
        [x1,y1] = ginput(1);
        
        % Now decide how to treat these values...
        switch optMode
            case 'sameSize'
                x1 = round(x1);
                y1 = round(y1);
        
                % Add in a second element a fixed distance away...
                xv = [x1 x1+fixVal];
                yv = [y1 y1+fixVal];
                
            case 'diffSize'
                [~,tmpX] = min(abs(dpn.fig.grid.opx - x1));
                [~,tmpY] = min(abs(dpn.fig.grid.opy - y1));
                
                xv = [x1 dpn.fig.grid.opx(tmpX + fixVal)];
                yv = [y1 dpn.fig.grid.opy(tmpY + fixVal)];        
        end
        
    otherwise
        % There is no otherwise        
end

% This is a way to quit implementation of the function
optAxLim = get(dpn.fig.ax.opt(1),'XLim');
if xv(1) < optAxLim(1) && xv(2) < optAxLim(1)
    set(f0,'Enable','on');
    return
end


% Now comes the hard part - we have to convert these regions / points into
% the equivalent pixels in the MS image
%[msX,msY,xv,yv] = convertOpt2MS(dpn,xv,yv);
[msX,msY,xv,yv] = xxxAnnotationConvertOpt2MS(dpn,xv,yv);

% Draw the things on the axes
[p,p2,p3] = patchDotAnnotation(dpn.fig,xv,yv,col,msX,msY);
   
% Anonymous function for coloured cell
colorgen = @(color,text) ['<html><table border=0 width=400 bgcolor=',...
    color,'><TR><TD>',text,'</TD></TR> </table></html>'];
col2 = colorgen(rgb2hex(col),int2str(usd));

% Place the important information into a cell
tmpAnno = {[p p2 p3] usd col col2 ['ID ' int2str(usd)] xv yv msX msY};

% Concatenate to the original table...
if isfield(dpn,'anno')
    dpn.anno = cat(1,dpn.anno,tmpAnno);
else
    dpn.anno = tmpAnno;
end

% Find only unique entries in terms of colour
unqCols = cell2mat(dpn.anno(:,2));
[~,b,~] = unique(unqCols);

% Update the table...
tabDat = dpn.anno(b,[4 5]);
set(man.tab,'Data',tabDat);

% Reenable all of the pushbuttons
set(f0,'Enable','on');

% Save to the guidata
guidata(fig.fig,dpn);

% Determine if the user is doing continuous annotation...
cont = get(man.cont,'Value');
if cont
    dpnDrawAnnotation(src,[],fig,man);
end
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,p2,p3] = patchDotAnnotation(fig,xv,yv,col,msX,msY)
% Place the patches / dots on the axes

% This bit actually draws the patch / dot for the annotated area
if xv(1) == xv(2) && yv(1) == yv(2)
    % Then we need to draw a single dot as the marker rather than a patched
    % rectangle...
    hold on;
    p = scatter(xv(1),yv(1),200,col,'o','filled',...
        'MarkerEdgeColor','w');
    
    % Now the MS axes...
    axes(fig.ax.ms1(1));
    hold on;
    p2 = scatter(msX,msY,200,col,'o','filled',...
        'MarkerEdgeColor','w');
    
    if isfield(fig.ax,'ms2')
        axes(fig.ax.ms2(1));
        hold on;
        p3 = scatter(msX,msY,200,col,'o','filled',...
            'MarkerEdgeColor','w');
    else
        p3 = [];
    end
else
    
    % Draw a patch on the optical axes...
    p = patch([xv(1) xv(2) xv(2) xv(1)],[yv(1) yv(1) yv(2) yv(2)],...
        col,...
        'EdgeColor',col,...
        'FaceColor',col,...
        'FaceAlpha',0.4,...
        'LineWidth',3);
    
    % Quicky need to modify msX and msY to extend themselves to the edge of
    % the pixels rather than only extending to their centres
    msX = msX + [-0.5 -0.5 0.5 0.5];
    msY = msY + [-0.5 0.5 0.5 -0.5];
    
    axes(fig.ax.ms1(1));
    hold on;
    p2 = patch(msX,msY,col,...
        'EdgeColor',col,...
        'FaceColor',col,...
        'FaceAlpha',0.4,...
        'LineWidth',3);
    
    % Only draw if in dual mode
    if isfield(fig.ax,'ms2')        
        axes(fig.ax.ms2(1));
        hold on;
        p3 = patch(msX,msY,col,...
            'EdgeColor',col,...
            'FaceColor',col,...
            'FaceAlpha',0.4,...
            'LineWidth',3);
    else
        p3 = [];
    end
        
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [px2,py2,x2,y2] = convertOpt2MS(dpn,x,y)
% Determine the MS location based on the region selected in the optical
% image

% First stage is to differentiate between a single click annotation and a
% rectangular annotation. This is required to ensure that the boundaries
% are specified to match the grid lines
if x(2) == x(1) && y(2) == y(1)
    type = 'single';
else
    type = 'rect';
end

% Now we need to convert the clicked optical pixels to match coordinates in
% the grid.  Note that the gridlines lie between pixels, rather than at the
% pixel centre (and its true coordinate)
switch type
    
    case 'rect'
        % Here we need to change the clicked pixels to match the grid, 
        % which will be the tricky part
        [~,x2(1)] = min(abs(dpn.fig.grid.opx - x(1)));
        [~,x2(2)] = min(abs(dpn.fig.grid.opx - x(2)));
        x2 = dpn.fig.grid.opx(x2);

        % Do the same for the y coordinates
        [~,y2(1)] = min(abs(dpn.fig.grid.opy - y(1)));
        [~,y2(2)] = min(abs(dpn.fig.grid.opy - y(2)));
        y2 = dpn.fig.grid.opy(y2);
        
        % Make sure that they are in the right order
        xx = sort(x2);
        yy = sort(y2);
        
    case 'single'
        % Find the nearest grid in the x and y dimension...
        [~,x2(1)] = min(abs(dpn.fig.grid.opx - x(1)));
        [~,y2(1)] = min(abs(dpn.fig.grid.opy - y(1)));
        
        % Now determine the other neighbour in the x dimension
        if dpn.fig.grid.opx(x2(1)) < x(1)
            x2(2) = x2(1) + 1;
        else
            x2(2) = x2(1) - 1;
        end
        
        % Same in the y dimension
        if dpn.fig.grid.opy(y2(1)) < y(1)
            y2(2) = y2(1) + 1;
        else
            y2(2) = y2(1) - 1;
        end
        
        % Now that we have the two grid lines then we take the mean to
        % ensure that the spot ends up in the middle...
        xx = mean(dpn.fig.grid.opx(x2));
        yy = mean(dpn.fig.grid.opy(y2));
        x2 = [xx xx];
        y2 = [yy yy];
        
end
        

% Shouldn't have to do this as it always finds the nearest gridline, rather
% than the clicked value.  This should bring wayward clicks back to be
% somewhere within the axes.  Perhaps...
% Ensure that xx and yy are within the optical boundaries...
%maskx = xx >= 1 & xx <= size(dpn.opt.coreg,2);
%masky = yy >= 1 & yy <= size(dpn.opt.coreg,1);
%xx = xx(maskx);
%yy = yy(masky);

% Notional coordinates for the optial image patch.  These should match the
% gridlines on the optical image
px = [xx(1) xx(1) xx(end) xx(end)];
py = [yy(1) yy(end) yy(end) yy(1)];

% Check for pixels outside the image boundary
px(px < 0) = 1;
py(py < 0) = 1;
px(px > size(dpn.opt.coreg,2)) = size(dpn.opt.coreg,2);
py(py > size(dpn.opt.coreg,1)) = size(dpn.opt.coreg,1);

% The resizing currently performed is not appropriate, as the effect is
% non-linear over the whole image.  Instead, the scaling factor should be
% proportional, as currently high numbered-pixels are trimmed
% disproportionally more than those with low pixel numbers
cx = px ./ size(dpn.opt.coreg,2);
cy = py ./ size(dpn.opt.coreg,1);

% The vales of cx and cy are relative ratios for annotation, so these just
% need to be translated to the MS image size
px2 = (cx .* size(dpn.d1.img,2));
py2 = (cy .* size(dpn.d1.img,1));

% Floor and ceil to maximise the annotated area...
if numel(xx) == 1 && numel(yy) == 1
    % Then we just need a single dot
    px2 = round(px2(1));
    py2 = round(py2(1));
else
    px2 = [ceil(px2(1)) ceil(px2(2)) floor(px2(3)) floor(px2(4))];
    py2 = [ceil(py2(1)) floor(py2(2)) floor(py2(3)) ceil(py2(4))];
end

% But have to check that values are within in the image boundary sizes
px2(px2 < 1) = 1;
py2(py2 < 1) = 1;
px2(px2 > size(dpn.d1.img,2)) = size(dpn.d1.img,2);
py2(py2 > size(dpn.d1.img,1)) = size(dpn.d1.img,1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%