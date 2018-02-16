function [px2,py2,x2,y2] = xxxAnnotationConvertOpt2MS(dpn,x,y)
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


