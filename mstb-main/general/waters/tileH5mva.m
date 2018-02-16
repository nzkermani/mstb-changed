function tileH5mva(gc,mva,comps)
%tileH5mva - tile the results of mva in an image, filing the empty gaps as
%required...

% What are we to plot in the images?
if nargin == 2
    comps = [1];
end

% Tiling arrangement.
numCol = 12;

% How many files are there?
numF = size(gc.fn,1);

% Determine image sizings
[numRow,pixRow,pixCol,sz] = detImgSize(numF,numCol,gc.sz);

% Now we can start to assemble the parts in each row
for r = 1:numRow
    
    % Create an image of the correct size...
    if comps(1) <= size(mva.ss,2)
        tmp = NaN(pixRow(r),pixCol,numel(comps));
    else
        tmp = NaN(pixRow(r),pixCol);
    end
    
    % Insert each part in to it...
    for c = 1:numCol
        
        % File index
        n = sz(r,c,1);
        if n == 0
            continue;
        end
        
        % Row/col index
        ri = [1 sz(r,c,2)];
        ci = [sum(sz(r,1:c,3))-sz(r,c,3)+1 sum(sz(r,1:c,3))+0];
        
        % Determine the image to be kept
        fx = gc.idx == n;
        
        if comps(1) <= size(mva.ss,2)
            
            % Select the correct pixels
            fy = mva.ss(fx,comps);
            
            % Set background pixels to NaN values
            tobg = gc.tobg(fx,1) == 0;
            fy(tobg,:) = NaN;
            
            % Final reshape into image format
            fy = reshape(fy,[sz(r,c,2) sz(r,c,3) numel(comps)]);
                        
        else
            % Draw an ion image instead
            mzf = mzFind(gc.mz,comps(1),comps(2));
            
            % This is the image
            fy = gc.sp(fx,mzf);
            
            % Set background pixels to NaN values
            tobg = gc.tobg(fx,1) == 0;
            fy(tobg,:) = NaN;
            
            % Reshape into an image
            fy = reshape(fy,[sz(r,c,2) sz(r,c,3)]);
                        
            %fy = bsxfun(@rdivide,fy,nansum(fy(:))) * 1000;
            
            % Normalise and log... although perhaps avoid here
            %fy = bsxfun(@rdivide,fy,max(fy(:))) * 1000;
            fy = log(fy + 1);
        end
        
        
        % Now dump in the image...
        tmp(ri(1):ri(2),ci(1):ci(2),:) = fy;
        
    end
    
    if r == 1
        img = tmp;
    else
        img = cat(1,img,tmp);
    end
end

% Here we might wish to change the image if we are doing kmeans
switch mva.method
    
    case 'kmeans'
        % Convert to colours instead...
        
        % Reshape
        sz = size(img);
        img = reshape(img,[sz(1)*sz(2) 1]);
        img(isnan(img)) = 0;
        
        % Unique entries
        [unq,~,ind] = unique(img(:));
        cols = parula(numel(unq)-1);
        
        % New image
        img2 = NaN(size(img,1),3);
        
        % Fill it up
        for n = 2:numel(unq)
            fx = ind == n;            
            img2(fx,:) = repmat(cols(n-1,:),[sum(fx) 1]);            
        end       
            
        % Restore original image size
        img = reshape(img2,[sz(1) sz(2) 3]);
            
        
    otherwise
        % Let's slightly modify the 0-5 and 95-100th percentiles, by altering their
        % intensities a little to make images a little easier on the eye...
        for n = 1:size(img,3)%numel(comps)       
            tmp = img(:,:,n);
            prc = prctile(tmp(:),[5 95]);        
            tmp(tmp < prc(1)) = prc(1);
            tmp(tmp > prc(2)) = prc(2);        
            img(:,:,n) = tmp;
        end
end

% If there are three components, then we are supposed to
% display this as a scaled image, hence we need to scale values
% between 0 and 1
if numel(comps) == 3
    img = imScale(img);
end

% Sort out the axes and stuff
figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);
imagesc(img);
axis image
axis off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [numRow,pixRow,pixCol,sz] = detImgSize(numF,numCol,imgSize)
% Determine image sizings

% Determine the arrangement of the various files, in rows of numCols
numRow = ceil(numF/numCol);

% Create a cell for sizes to determine the array size...
sz = zeros(numRow,numCol,3);
for n = 1:numF
    
    % Determine row and col placements for this file...
    row = ceil(n/numCol);
    col = mod(n,numCol);
    if col == 0
        col = numCol;
    end
    
    % Place the index in the sz matrix (1)
    sz(row,col,1) = n;
    
    % Place the sizes in it also (2,3)
    sz(row,col,2) = imgSize(n,1);
    sz(row,col,3) = imgSize(n,2);
    
end
    
% Determine the maximum width of each tiled row
pixRow = max(sz(:,:,2),[],2);
pixCol = max(sum(sz(:,:,3),2));


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%