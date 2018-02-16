function [sp,method,isInterp] = dpnPerformInterpolation(sp,method,fx)
%dpnPerformInterpolation - from the split image we need to determine
%intensities in the pixels that were not measured

wb = waitbar(0,'Interpolating');

if nargin == 1
    method = 'linear';
end

% Ensure that a method has been selected...
switch method
    case {'linear','cubic','spline','nearest'}
        % Do nothing
        interpFunction = 'interp2';
        
    case {'bilinear','bicubic'}
        % Then we have an entirely different method of interpolation
        interpFunction = 'imresize';
        
    otherwise
        % Set to linear by default
        method = 'linear';
        interpFunction = 'interp2';
end

% Determine the full / empty rows in this instance. fx is a logical vector
% denoting the ORIGINALLY measured pixels/rows. Initially it isn't provided
% as it is worked out here. When performing interpolation again, then this
% information must be provided to ensure that we re-interpolate the
% original pixels, rather than using the interpolated data
if isempty(fx)
    totImg = nansum(sp,3);
    rowSum = nansum(totImg,2);

    % Find the full and the empty rows
    fx = rowSum > 0;
else
    % Totally provided, so need to do nothing
end

% This is a simple vector from 1 to n
yy = 1:size(sp,1);

% Orig coordinates
ox = 1:1:size(sp,2);
oy = yy(fx);
[mx,my] = meshgrid(ox,oy);

% New coordinates
%nx = mx;
%ny = yy;
[nx,ny] = meshgrid(ox,yy);

% Now we do some interpolation... start with totImg to see if it looks good
% as a starting point
numI = size(sp,3);

% Interpolate over each image
switch interpFunction
    case 'interp2'
        for n = 1:numI    
            sp(:,:,n) = interp2(mx,my,sp(fx,:,n),nx,ny,method);

            wb = waitbar(n/numI,wb,'Interpolating');
        end

    case 'imresize'
        for n = 1:numI    
            sp(:,:,n) = imresize(sp(fx,:,n),[size(sp,1) size(sp,2)],method);

            wb = waitbar(n/numI,wb,'Interpolating');
        end

        
end

% Create a matrix of pixel identity that says whether the pixels were
% original intensities (i.e. measured) or actually interpolated
isInterp = repmat(~fx,[1 size(sp,2)]);

% Finally close the waitbar
delete(wb);

end

