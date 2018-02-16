function [ img ] = desiIonTissueTile(im,imTIC,imMask,comps)
%desiIonTissueTile - visualisation of the images extracted from the other
%function: desiIonTissueExtract.m...

% Tiling arrangement.
numCol = 10;

% How many files are there?
numF = size(im,1);

% What are the sizes of the things?
[szImg] = sizeOfInputs(im);

% Determine image sizings
[numRow,pixRow,pixCol,sz] = detImgSize(numF,numCol,szImg);

% Now we can start to assemble the parts in each row
for r = 1:numRow
    
    % Create an image of the correct size...
    if numel(comps) == 2
        tmp = NaN(pixRow(r),pixCol,1);
    else
        tmp = NaN(pixRow(r),pixCol,numel(comps));
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

        % This is the image
        fy = im{n,1}(:,:,comps);
        
        % Normalise the intensities according to the TIC image...
        fy = bsxfun(@rdivide,fy,imTIC{n,1});
        
        
        % If there are two components, then we divide first by the second
        % and don't alter it...
        if size(fy,3) == 2
            [fy] = image2comps(fy);
        elseif size(fy,3) == 3
            fy = imScale(fy);
        end
        
        % Suppress low/high intensities for each channel of the image
        [fy] = intensitySuppress(fy,5,95);
        
        % Remove pixels which weren't those of the annotated tissue
        fy = bsxfun(@times,fy,imMask{n,1});
        
        % Now dump in the image...
        tmp(ri(1):ri(2),ci(1):ci(2),:) = fy;
        
    end
    
    % Combine into a final image...
    if r == 1
        img = tmp;
    else
        img = cat(1,img,tmp);
    end
end

% Draw the image
figure; imagesc(img);
axis image;



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [szImg] = sizeOfInputs(im)
% Determine the sizes of the individual data sets

szImg = cellfun(@size,im(:,1),'UniformOutput',false);

% Reshape into rows and colums rather than a cell array
if numel(szImg{1}) == 2
    szImg = reshape([szImg{:}],[2 size(im,1)])';
else
    szImg = reshape([szImg{:}],[3 size(im,1)])';
end

end
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
function [img] = image2comps(img)

img = (img(:,:,1) ./ img(:,:,2));
img(img == 0) = NaN;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [img] = intensitySuppress(img,prcLo,prcHi)

% Size of the image
sz = size(img,3);

for x = 1:sz
    
    % Temporarily extract image
    pctmp = img(:,:,x);
    
    % Determine the percentiles
    prc = prctile(pctmp(:),[prcLo prcHi]);
    
    % Error handling in case of zero values
    if sum(prc) == 0
        prc = prctile(pctmp(pctmp > 0),[5 95]);
    end
    
    % Set values outside of prcLo--prcHi to corresponding values
    pctmp(pctmp < prc(1)) = prc(1);
    pctmp(pctmp > prc(2)) = prc(2);
    
    % Return channel to the image
    img(:,:,x) = pctmp;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
