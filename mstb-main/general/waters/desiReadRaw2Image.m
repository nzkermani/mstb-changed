function [ data ] = desiReadRaw2Image(sp,xy2D)
%desiReadRaw2Image - put the cell array into an image-like cell array

% Determine the xy coordinates for the sample
[nRows,nCols] = size(xy2D);

% Preallocate the spectral datacube. Can we make this sparse without 
% affecting the rest of the code?
data = cell(nRows,nCols);
%tt = zeros(nRows,nCols);

% Loop through each scan and align the m/z vectors
for x = 1:nRows
    for y = 1:nCols
        if ~isnan(xy2D(x,y))
            
            % Coords
            idx = xy2D(x,y);
            
            % Info
            tmp = sp{idx,1}';
            
            data{nRows-x+1,y} = tmp;
            
            %tt(nRows-x+1,y) = sum(tmp(:,2));
        end
            
    end
end

%figure; imagesc(tt);

end

