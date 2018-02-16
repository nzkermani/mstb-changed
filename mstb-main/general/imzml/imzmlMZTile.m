function [ img ] = imzmlMZTile(im,comps,mz)
%imzmlMZTile - tile the imzML images from the imzmlIonExtractBatch
%function.

% Tiling arrangement.
numCol = 15;

% How many files are there?
numF = size(im,1);

% What are the sizes of the things?
szImg = cellfun(@size,im(:,1),'UniformOutput',false);
if numel(szImg{1}) == 2
    szImg = reshape([szImg{:}],[2 size(im,1)])';
else
    szImg = reshape([szImg{:}],[3 size(im,1)])';
end

% Determine image sizings
[numRow,pixRow,pixCol,sz] = detImgSize(numF,numCol,szImg);

% Now we can start to assemble the parts in each row
for r = 1:numRow
    
    % Create an image of the correct size...
    tmp = NaN(pixRow(r),pixCol,numel(comps));
    
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
        %fx = gc.idx == n;
        
        % Draw an ion image instead
        %mzf = mzFind(gc.mz,comps(1),comps(2));

        % This is the image
        fy = im{n}(:,:,comps);

        % Normalise and log... although perhaps avoid here
        %fy = bsxfun(@rdivide,fy,max(fy(:))) * 1000;
        %fy = log(fy + 1);
        
        % Trim the percentiles to make the image stand out a little better
%         for x = 1:size(fy,3)%numel(comps)       
%             pctmp = fy(:,:,x);
%             prc = prctile(pctmp(:),[5 95]);        
%             pctmp(pctmp < prc(1)) = prc(1);
%             pctmp(pctmp > prc(2)) = prc(2);        
%             fy(:,:,x) = pctmp;
%         end
%         
%         % Scale the image parts between 0-1
%         fy = imScale(fy);
        
        %fy(isnan(fy)) = 0;
        
        % Now dump in the image...
        tmp(ri(1):ri(2),ci(1):ci(2),:) = fy;
        
    end
    
    if r == 1
        img = tmp;
    else
        img = cat(1,img,tmp);
    end
end

% Divide by the expected m/z value
img = 1e6 * (img - mz) / mz;

figure; imagesc(img);
axis image;
colorbar
colormap('redbluecmap');

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

