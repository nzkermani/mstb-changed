function [mz,data] = desiRawPeakExtract(pks,img)
%desiRawPeakExtract - use the peaks / wdiths to extract images

% Sizes!
numP = numel(pks.mz);
sz = size(img);

% Create empty image
data = zeros(sz(1),sz(2),numP);
mzs = NaN(size(data));

wb = waitbar(0);

% Loop through each pixel in turn...
for x = 1:sz(1)
    
    for y = 1:sz(2)
        
        if isempty(img{x,y})
            continue;
        end
        
        % Now loop through each peak
        for n = 1:numP
            
            mzLo = pks.mz(n) - pks.hw(n);
            mzHi = pks.mz(n) + pks.hw(n);
            
            % Mask
            fx = img{x,y}(:,1) >= mzLo & img{x,y}(:,1) <= mzHi;
            
            % What if we find multiple ones?
            if sum(fx) > 1
                % Should we use the one closest to the centroid?
                tmp1 = img{x,y}(fx,1);
                tmp2 = img{x,y}(fx,2);                
                [~,b] = min(abs(tmp1 - pks.mz(n)));                
                data(x,y,n) = tmp2(b);
                
                mzs(x,y,n) = tmp1(b);
                
            elseif sum(fx) == 1
                % Save into matrix...
                data(x,y,n) = img{x,y}(fx,2);
                mzs(x,y,n)  = img{x,y}(fx,1);
                
            end
                        
        end
        
    end
    
    waitbar(x/sz(1),wb);
    
end
           
delete(wb);

% Format the mz values, which are calculated from the median of the raw
% data, rather than the local maxima or anything else like that
mz = nanmedian(reshape(mzs,[sz(1)*sz(2) numP]),1);

            


end

