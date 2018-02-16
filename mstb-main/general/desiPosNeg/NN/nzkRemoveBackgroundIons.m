function [ chk ] = nzkRemoveBackgroundIons(sp,tobg,tobgImg)
%nzkRemoveBackgroundIons - just determine the ions which are in the
%background

% Determine the mean intensities over true|false values
mean0 = mean(tobgImg(tobg == 0));
mean1 = mean(tobgImg(tobg == 1));

% Only proceed if mean1 >> mean0
if mean1 > (mean0 * 2)
    
    chk = false(size(sp,3),1);
    for n = 1:size(sp,3)
        
        tmp = squeeze(sp(:,:,n));
        m0 = nanmean(tmp(tobg == 0));
        m1 = nanmean(tmp(tobg == 1));
        
        if m1 > (m0 * 2)
            chk(n,1) = true;
        end
    end
    
else
    % Cannot be sure that the tobg segmentation is correct
    disp('Cannot remove background ions');
    chk = true(size(sp,3),1);
    
end

end

