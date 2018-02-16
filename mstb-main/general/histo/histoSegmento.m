function [idx] = histoSegmento(img)
%histoSegmento - does a little bit of clustering on H&E images. Want to see
%if it is possible to do histology-driven annotation. Certainly start with
%differentiation of tissue v non-tissue object pixels. Use the high-res
%optical image to make it more impressive. Can always be done low-res if
%required.
%
% James McKenzie, 2016
%
% Reference
% Much of the background to H&E segmentation is taken from here:
% http://uk.mathworks.com/help/images/examples/
%   color-based-segmentation-using-k-means-clustering.html
%

% Some default options
method = 'kmeans';
maxClust = 3;
doTOBG = false;

% May wish to smooth the image a little - which is probably a lot better
% than reducing the resolution
%
%

% Convert to a new colour space - see link below
[img] = imgConvert(img);

% Determine TO/BG
if doTOBG
    [tobg] = imgTOBG(img);
else
    tobg = ones(size(img,1)*size(img,2),1);
end
    

% Run the clustering function
[idx] = performClustering(img,tobg,method,maxClust);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [img] = imgConvert(img)
% Change the colourspace of the image from RGB to CIE L*a*b

cform = makecform('srgb2lab');
img = applycform(img,cform);

% Convert to a double?
img = double(img);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tobg] = imgTOBG(img)

% Scale
img = img ./ max(img(:));

% Convert this to grayscale
gray = nansum(img,3);
gray = max(gray(:)) - gray;

% Run the Otsu tresholding method
[tobg,~] = dpnTOBG(gray,[],[]);

% Perhaps we can consider smoothing the image here...
filt = fspecial('average',10);
tobg = filter2(filt,tobg) >  0.5;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [idx] = performClustering(img,tobg,method,maxClust)

% Size of the image
sz = size(img);

% Reshape the image
img = reshape(img, [sz(1)*sz(2) sz(3)]);

% And the tobg map with a logical vector
tobg = reshape(tobg,[sz(1)*sz(2) 1]) == 1;

% We only need the 'a' and 'b' parts of the colour scheme (and TO not BG)
img = img(tobg,[2 3]);

% Here we do one of some clustering methods
idx = zeros([sz(1) sz(2) maxClust-1]);
i = 0;
for n = 2:maxClust
    
    % Counter
    i = i + 1;
    
    % Perform clustering
    switch method
        
        case {'k','kmeans'}
              
            % k-means clustering
            [tmp,cent] = kmeans(img(tobg,:),n,...
                'Distance','sqEuclidean',...
                'Replicates',1);
            
            
        case {'g','gauss','gaussian'}
            
            % Gaussian mixture models
            options = statset('Display','final');
            gm = fitgmdist(img,n,'Options',options);
            
            
            % ri = round(linspace(1,size(img,1),1000));
            % figure; hold on;
            % scatter(img(ri,1),img(ri,2));
            % ezcontour(@(x,y)pdf(gm,[x y]),xlim,ylim);
            % drawnow;
            
            % Need to sort out the clusters...
            tmp = cluster(gm,img);
            
            pro = posterior(gm,img);
            % fx = sum(tmp,2) > 1;
            % tmp(fx,:) = false;
            % tmp = bsxfun(@times,tmp,[1:size(tmp,2)]);
            % tmp = sum(tmp,2);
            probs = zeros(sz(1)*sz(2),size(pro,2));
            probs(tobg,:) = pro;
            probs = reshape(probs,[sz(1) sz(2) size(pro,2)]);
            
        otherwise
            % There is no otherwise
    end
    
    % Place back in original locations considering the tobg pixels
    tmp2 = zeros(sz(1)*sz(2),1);
    tmp2(tobg,1) = tmp;
    
    % Reshape and put in the results matrix
    idx(:,:,i) = reshape(tmp2,[sz(1) sz(2)]);
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%