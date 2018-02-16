function [d1,d2,anno] = dpnPrepFullMS(dpn,mzRange,tobg)
%dpnPrepFullMS - get the data for each MS mode and prepare it in an
%appropriate matrix.

% Trim here
[d1.sp,d1.idx,d1.mz] = getDXmatrix(dpn.d1.sp,dpn.d1.mz,mzRange,tobg);
d1.mmcClass = [];
d1.mmcMask2 = [];

% Set all NaN values to 0. Not sure why these have persisted if I'm honest.
% This is likely to be an interpolation issue.
mask = isnan(d1.sp);
d1.sp(mask) = 0;

% Do the same in dual mode if necessary
if isfield(dpn,'d2')
    [d2.sp,d2.idx,d2.mz] = getDXmatrix(dpn.d2.sp,dpn.d2.mz,mzRange,tobg);
    mask = isnan(d2.sp);
    d2.sp(mask) = 0;
    d2.mmcClass = [];
    d2.mmcMask2 = [];

else
    d2 = [];
end

% Here we extract the annotations.
if strcmp(dpn.mode,'dual')
    [anno.mask2,...
        anno.histID,...
        anno.isInt1,...
        anno.isInt2,...
        anno.pixID] = dpnAnnotationExtract(dpn);

elseif strcmp(dpn.mode,'single')
    [anno.mask2,...
        anno.histID,...
        anno.pixID] = desiAnnotationExtract(dpn);
end
    
if ~isempty(tobg)
    tobg = reshape(tobg,[numel(tobg) 1]);
    
    anno.mask2 = anno.mask2(tobg);
    anno.histID = anno.histID(tobg);
    anno.pixID = anno.pixID(tobg);
    
    if strcmp(dpn.mode,'dual')
        anno.isInt1 = anno.isInt1(tobg);
        anno.isInt2 = anno.isInt2(tobg);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx,idx,newMZ] = getDXmatrix(sp,mz,mzRange,tobg)
% From a 3D image, strip out background if requested, and return to a 2D
% matrix with which to do PCA and stuff...

% Determine the m/z mask that defines which variables to use
mask = mz > min(mzRange) & mz < max(mzRange);

newMZ = mz(mask);

% Size of sp
sz = size(sp);

% Do we need to reshape tobg?
if ~isempty(tobg)    
    flag = true;
    tobg = reshape(tobg,[sz(1)*sz(2) 1]);    
else
    flag = false;
end

% Reshape the data matrix
dx = reshape(sp,[sz(1)*sz(2) sz(3)]);
dx = dx(:,mask);

% This is necessary for use when resizing background removed images
idx = [1:size(dx,1)]'; %#ok<NBRAK>

% Remove tobg?
if ~flag
    return
end

% Trim out the background here
dx = dx(tobg == 1,:);
idx = idx(tobg == 1,:);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%