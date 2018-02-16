function [mask2,histID,isInt1,isInt2,pixID] = dpnAnnotationExtract(dpn)
%dpnAnnotationExtract - need to extract indices for each data set

% For the MS images, we need to create a blank mask the same size as the
% image.  It is zero by default, and non-zero if the pixels have been
% annotated.  We will place the usd value into the matrix for each pixel.
% Once this has been accomplished then we can reshape this matrix and
% both of the MS images (pos/neg). Then pixels are extracted and a histID
% classification vector is determined. Then we can run the stats function
% or something similar
sz = size(dpn.d1.img);
numA = size(dpn.anno,1);

mask = zeros(sz);
histID = cell(sz);
isInt1 = cell(sz);
isInt2 = cell(sz);
pixID  = cell(sz);

% Now loop through each of the annotations
for n = 1:numA
    
    x = dpn.anno{n,8};
    y = dpn.anno{n,9};
    
    for i = min(x):1:max(x)
        for j = min(y):1:max(y)
            mask(j,i) = dpn.anno{n,2};
            histID(j,i) = dpn.anno(n,5); 
            pixID(j,i) = {['Row/Col ' int2str(j) '-' int2str(i)]};
            
            % Interpolation performed?
            tmp = dpn.d1.isInterp(j,i);
            if tmp
                isInt1(j,i) = {'Interpolated'};
            else
                isInt1(j,i) = {'Measured'};
            end
            
            % For negative too
            tmp = dpn.d2.isInterp(j,i);
            if tmp
                isInt2(j,i) = {'Interpolated'};
            else
                isInt2(j,i) = {'Measured'};
            end
        end
    end
    
end

mask2 = reshape(mask,[sz(1)*sz(2) 1]);
histID = reshape(histID,[sz(1)*sz(2) 1]);
isInt1 = reshape(isInt1,[sz(1)*sz(2) 1]);
isInt2 = reshape(isInt2,[sz(1)*sz(2) 1]);
pixID  = reshape(pixID, [sz(1)*sz(2) 1]);

end

