function [ tobg ] = pixTOBG(img,nBins,smooth)
%pixTOBG - determine TO / BG pixels, and then possibly smooth the thing too
%in order to get a rounded image

% This function runs the Otsu thresholding approach
if isempty(nBins)
    [tobg,nBins] = getObjectPixels(img);
else
    [tobg,nBins] = getObjectPixels(img,nBins);
end

% Now maybe smooth the image...
if ~isempty(smooth)
    
    % First smooth to determine the boundaries of the tissue section
    filt = fspecial('average',3);
    tobg = filter2(filt,tobg) >  0.75;
    
    % Second smooth to close off insides
    %smo2 = filter2(filt,smo1);
    
end


%figure; 
%ax1 = subplot(1,3,1); imagesc(tobg);
%ax2 = subplot(1,3,2); imagesc(smo1);
%ax3 = subplot(1,3,3); imagesc(smo2);



    


end

