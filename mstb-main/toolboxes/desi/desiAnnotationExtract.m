function [mask2,histID,pixID] = desiAnnotationExtract(dpn,flag)
%desiAnnotationExtract - need to extract indices for each data set

% For the MS images, we need to create a blank mask the same size as the
% image.  It is zero by default, and non-zero if the pixels have been
% annotated.  We will place the usd value into the matrix for each pixel.
% Once this has been accomplished then we can reshape this matrix and
% the MS image. Then pixels are extracted and a histID
% classification vector is determined. Then we can run the stats function
% or something similar
sz = [size(dpn.d1.sp,1) size(dpn.d1.sp,2)];
numA = size(dpn.anno,1);

mask = zeros(sz);
histID = cell(sz);
pixID  = cell(sz);
annoID = zeros(sz);

% Now loop through each of the annotations
for n = 1:numA
    
    x = dpn.anno{n,8};
    y = dpn.anno{n,9};
    
    for i = min(x):1:max(x)
        for j = min(y):1:max(y)
            mask(j,i) = dpn.anno{n,2};
            histID(j,i) = dpn.anno(n,5); 
            pixID(j,i) = {['Row/Col ' int2str(j) '-' int2str(i)]};
            annoID(j,i) = n;
        end
    end
    
end

mask2 = reshape(mask,[sz(1)*sz(2) 1]);
histID = reshape(histID,[sz(1)*sz(2) 1]);
pixID  = reshape(pixID, [sz(1)*sz(2) 1]);
annoID = reshape(annoID,[sz(1)*sz(2) 1]);

% Output the annotation ID if specified
if nargin == 1
    flag = false;
end
if ~islogical(flag)
    flag = false;
end
if flag
    pixID = annoID;
end

end

