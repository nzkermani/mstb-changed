function dpnSegmentationPerform(~,~,fig,seg)
%dpnSegmentationPerform - calculate the segmentatation requested and work
%out how to show it in the figures...

tic

% Get the data
dpn = guidata(fig.fig);

% Get all of the input parameters
[~,~,~,mzL,mzH,choice,tobg,fuse,maxClust] = getAllOptions(dpn,seg);

% Determine if we can continue - there is an option (perhaps unnecessary)
% to quit
if strcmpi(choice,'none')
    redrawIonImages(dpn);
    return
end

% Extract the data matrices for both images.  This currently works for the
% PCA/kMeans/MMC analyses, though more will be added in the future.
switch dpn.mode
    
    case 'dual'
        dx(3,1) = struct('sp',[],'idx',[],'mz',[],...
            'mmcClass',[],'mmcMask2',[]);
        [dx(1,1),dx(2,1),anno] = dpnPrepFullMS(dpn,[mzL mzH],tobg);
        
    case 'single'
        [dx,~,anno] = dpnPrepFullMS(dpn,[mzL mzH],tobg);      
end
      
% Create a waitbar.
wb = waitbar(0,'Performing...');

% Let's decide first about the fusion method - if high level then we need
% to do norm/etc on each dataset prior to a PCA. If low level, then we
% concatenate and normalise together...
if strcmp(dpn.mode,'dual')
    [dx] = prepareFusion(dx,seg,fuse);
end

% Here is where we can be specific about the analysis to be performed. We
% have extracted the variables to be used above, so here is where we
% diversify
switch choice
    case 'PCA'
        [dx] = segMethodPCA(dx,wb);
        
    case 'kMeans'
        [dx] = segMethodKmeans(dx,wb,maxClust);
        
    case 'MMC'
        [dx] = xxxSegMethodMMC(dx,wb,anno);
        
    otherwise
        % There is no otherwise
        disp('There is no otherwise');
        return
end

% Add the annotated labels
%dx.anno = unique(anno.histID(anno.mask2 > 0));

% Return the 'data' element of the structures to images!
[dpn] = generalPreparationFunction(dpn,dx,choice);

% We need to ensure that the images are consistent in the ± respect. We
% might need to invert some of the PCs to make them visually compatible
if strcmp(dpn.mode,'dual')
    [dpn] = ensureImageCompatibility(dpn,choice);
end

% Save to the guidata
guidata(fig.fig,dpn);

% Delete the waitbar
delete(wb);

% From here it is all about presentation of the results!
resultsDisplayOptions(dpn,seg,choice);

toc

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [norm,doLog,tran,mzL,mzH,...
    choice,tobg,fuse,maxClust] = getAllOptions(dpn,seg)
% Determine all of the input parameters from the window

% Normalisation
methods = get(seg.norm,'String');
val = get(seg.norm,'Value');
norm = methods{val};

% Do we apply a log transform to the data?
doLog = get(seg.log,'Value') == 2;

% Scaling
methods = get(seg.scale,'String');
val = get(seg.scale,'Value');
tran = methods{val};

% Fusion method
if strcmp(dpn.mode,'dual')
    fv = get(seg.fuse,'Value');
    ft = get(seg.fuse,'String');
    fuse.method = ft{fv};
    if strcmpi(fuse.method(1:2),'hl')
        fuse.value = str2double(get(seg.fuseComp,'String'));
        fuse.retroLoad = true(get(seg.retroLoadings,'Value'));
        fuse.retroScores = true(get(seg.retroScores,'Value'));
    else
        fuse.value = [];
        fuse.retroLoad = [];
        fuse.retroScores = [];
    end
else    
    % Single mode has no fusion
    fuse.value = [];
    fuse.retroLoad = [];
    fuse.retroScores = [];
end

% Also determine the m/z values to be used
mzL = str2double(get(seg.mzL,'String'));
mzH = str2double(get(seg.mzH,'String'));

% Now determine which of the methods is to be performed
methods = get(seg.list,'String');
val = get(seg.list,'Value');
choice = methods{val};

% Maximum number of clusters
try
    maxClust = str2num(get(seg.maxClust,'String')); %#ok<ST2NM>
catch
    maxClust = 3;
    set(seg.maxClust,'String','3');
end

% Ignore the background pixels?
bg = get(seg.remBG,'Value');
if bg
    if isfield(dpn.d1,'tobg')
        tobg = dpn.d1.tobg;
    elseif isfield(dpn,'d2')
        if isfield(dpn.d2,'tobg')
            tobg = dpn.d2.tobg;
        else
            tobg = [];
        end
    else
        tobg = [];
    end    
    
    % If there is nothing, then just determine it automatically and be done
    % with it
    if isempty(tobg)
        [dpn.d1.tobg,~] = dpnTOBG(dpn.d1.img,[],[]);
        tobg = dpn.d1.tobg;
    end
else
    tobg = [];
end

% We need to fill in the background image quite a lot in order to make it
% one homogeneous area rather than potentially having a lot of holes in it
if ~isempty(tobg)
    
    % Just define the options here - may need to make these interactive
    % when we move towards different sections
    %mOpts.imopen = 2;
    %mOpts.imclose = 8;
    %mOpts.bigoper = [];
    
    % This does the operations defined.
    %tobg = doMorph(tobg,mOpts);
    tobg = tobg == 1;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx] = prepareFusion(dx,seg,fuse)
% Do the preparation for the fusion method.

% Decide what to do...
switch lower(fuse.method(1:2))
    
    case 'll'
        
        % Combine the two matrices together
        dx(3).sp = [dx(1).sp dx(2).sp];
        dx(3).idx = dx(1).idx;
        dx(3).mz = [dx(1).mz dx(2).mz];
        
        % Now we can do the norm/etc on the three matrices
        for n = 1:3
            [dx(n).mz,dx(n).sp] = xxxNormTranScal(seg,dx(n).mz,dx(n).sp);
        end
        
    case 'hl'
                
        % Here we do norm/etc on the two images
        for n = 1:2
            [dx(n).mz,dx(n).sp] = xxxNormTranScal(seg,dx(n).mz,dx(n).sp);
        end
        
        % Now we run PCA on each of these and select the number of
        % components dependent on the option
        [dx(3).mz,dx(3).sp,~] = xxxFusionHighLevel(dx(1:2),...
            fuse.method,...
            fuse.value);
        
        
    otherwise
        disp('There is no otherwise');
        return
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx] = segMethodPCA(dx,wb)
% Run PCA on the data sets

numD = max(size(dx));

for n = 1:numD
    
    % Generic waitbar progress
    waitbar(n/numD,wb,'PCA');

    % Calculate PCs
    [dx(n).ll,dx(n).data,dx(n).ee] = pca(dx(n).sp,...
        'NumComponents',10);    
    
    % Convert eigenvalues
    dx(n).ee = 100 * dx(n).ee(1:10) / sum(dx(n).ee);
    
    % Create labels...
    dx(n).labs = cell(10,1);
    for r = 1:10
        dx(n).labs(r,1) = {['PC' int2str(r) ', ' sprintf('%0.2f',dx(n).ee(r)) '%']};
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx] = segMethodKmeans(dx,wb,maxK)
% Run the kmeans clustering technique

% How many datasets?
numD = max(size(dx));

for n = 1:numD
    
    % Generate empty matrices for results
    sz = size(dx(n).sp,1);
    dx(n).data = zeros(sz,maxK-1);
    
    % Empty labels...
    dx(n).labs = cell(maxK-1,1);
    
    for k = 2:maxK
        
        % Waitbar
        waitbar((k-1)/(maxK-1),wb,[int2str(n) '/' int2str(numD) ': k = ' int2str(k)]);
        
        dx(n).data(:,k-1) = kmeans(dx(n).sp,k);
        
        % Write the little label
        dx(n).labs(k-1,1) = {['k = ' int2str(k)]};
        
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx] = segMethodMMC(dx,wb,anno)
% Here is the MMC method

% How many datasets?
numD = max(size(dx));

% Indices of the training set
idx = anno.mask2 > 0;

% Loop through them
for n = 1:numD    
    
    % What is the mean spectrum?
    meanSpec = nanmean(dx(n).sp(idx,:),1);

    % Subtract the mean.
    tmp1 = bsxfun(@minus,dx(n).sp(idx,:),meanSpec);
    
    % Run the OAA-MMC function
    [B1,W1] = oaaMMC(tmp1,anno.histID(idx),2);
    
    % Calculate probabilities for pixels
    [~,Pr1,~] = calcScoresSpectra(dx(n).sp,meanSpec,W1,B1);
    
    % Turn unambiguous probabilities into a vector of histological
    % classification (not used in this function)
    probs = Pr1 > 0.95;
        
    % Set ambiguous ones to zero
    check = sum(probs,2) ~= 1;
    probs(check,:) = 0;

    % Just give me a vector of class membership
    class = bsxfun(@times,double(probs),1:size(probs,2));
    class = max(class,[],2);
    
    % Cell of annotations
    mmcClass = cell(size(probs,1),1);
    mmcNumbs = zeros(size(probs,1),1);
    mmcAnnos = unique(anno.histID(anno.mask2 > 0));
    for r = 1:numel(mmcAnnos)        
        fx = class == r;
        mmcClass(fx,1) = mmcAnnos(r);
        
        % Add the annotated number too
        nu = max(anno.mask2(fx));
        mmcNumbs(fx,1) = nu;
    end

    % Create some labels
    dx(n).labs = cell(size(Pr1,2)+1,1);
    for r = 1:size(Pr1,2)
        dx(n).labs(r,1) = {['LV ' int2str(r)]};
    end
    dx(n).labs(r+1,1) = {'Composite'};
    
    % Save the probabilities to the structure
    dx(n).data = Pr1;
    
    % Save these two additional things that we have made
    dx(n).mmcClass = mmcClass;
    dx(n).mmcMask2 = mmcNumbs;
    
    % Update the waitbar
    waitbar(n/numD,wb,['MMC: ' int2str(n) '/' int2str(numD)]);
    
    % Information about this classification
    hh = hist(mmcNumbs,unique(mmcNumbs))
    sum(hh)
    100 * hh / sum(hh)

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function redrawIonImages(dpn)
% Do this when the user selects 'none'

% Return the images to the ion images...
set(dpn.fig.ax.ms1(2),'CData',dpn.d1.img);
set(dpn.fig.ax.ms2(2),'CData',dpn.d2.img);
set(dpn.fig.ax.ms1(1),...
    'XLim',[1 size(dpn.d1.img,2)],...
    'YLim',[1 size(dpn.d1.img,1)],...
    'XTick',[],...
    'YTick',[]);
set(dpn.fig.ax.ms2(1),...
    'XLim',[1 size(dpn.d1.img,2)],...
    'YLim',[1 size(dpn.d1.img,1)],...
    'XTick',[],...
    'YTick',[]);

dpnUpdateFusedIonImage([],[],dpn.fig);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dpn] = generalPreparationFunction(dpn,dx,choice)
% Return the data to images and then back into the dpn structure

mmcFlag = any(strfind(choice,'mmc')) || any(strfind(choice,'MMC'))

 sz  = size(dpn.d1.sp);
idx = dx(1).idx;
for n = 1:max(size(dx))
    
    % Why do we need to force 'real'?
    tmp = real(dx(n).data);
    
    % Create the temporary image
    tmpImg = xxxReturn2Image(tmp,idx,sz);
    
    % Determine what to do with it?
    if n == 1
        dpn.d1.mva.(choice).img = tmpImg;
        dpn.d1.mva.(choice).label = dx(n).labs;        
        if mmcFlag
            dpn.d1.mva.(choice).mmcClass = dx(n).mmcClass;
            dpn.d1.mva.(choice).mmcMask2 = dx(n).mmcMask2;
        end
        
    elseif n == 2
        dpn.d2.mva.(choice).img = tmpImg;
        dpn.d2.mva.(choice).label = dx(n).labs;        
        if mmcFlag
            dpn.d2.mva.(choice).mmcClass = dx(n).mmcClass;
            dpn.d2.mva.(choice).mmcMask2 = dx(n).mmcMask2;
        end
        
    elseif n == 3
        dpn.fuse.mva.(choice).img = tmpImg;
        dpn.fuse.mva.(choice).label = dx(n).labs;
        if mmcFlag
            dpn.fuse.mva.(choice).mmcClass = dx(n).mmcClass;
            dpn.fuse.mva.(choice).mmcMask2 = dx(n).mmcMask2;
        end
        
    else
        warning('Nothing else');
        return
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tmpData] = return2image(data,index,sz)

% What is the size of the new images?
newSz = [sz(1) sz(2) size(data,2)];

% Full size empty matrix
tmpData = NaN([newSz(1)*newSz(2) newSz(3)]);

% Restore to original location
tmpData(index,:) = data;

% Reshape
tmpData = reshape(tmpData,newSz);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dpn] = ensureImageCompatibility(dpn,choice)

switch choice
    case 'PCA'
        [dpn.d1,dpn.d2] = pcaImageCheck(dpn.d1,dpn.d2,choice);
            
    case 'kMeans'
        dpn.d1.mva.(choice).img(isnan(dpn.d1.mva.(choice).img)) = 0;
        dpn.d2.mva.(choice).img(isnan(dpn.d2.mva.(choice).img)) = 0;    
        [dpn.d1,dpn.d2] = kMeansImageCheck(dpn.d1,dpn.d2,choice);
        [dpn.d1,dpn.fuse] = kMeansImageCheck(dpn.d1,dpn.fuse,choice);
        
    case 'MMC'
        dpn.d1.mva.(choice).img(isnan(dpn.d1.mva.(choice).img)) = 0;
        dpn.d2.mva.(choice).img(isnan(dpn.d2.mva.(choice).img)) = 0;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx1,dx2] = pcaImageCheck(dx1,dx2,choice)

numI = size(dx1.mva.(choice).img,3);

for n = 1:numI
    
    % Check addition and subtraction of the images
    a1 = dx1.mva.(choice).img(:,:,n) + dx2.mva.(choice).img(:,:,n);

    s1 = dx1.mva.(choice).img(:,:,n) - dx2.mva.(choice).img(:,:,n);
    
    hst = hist([a1(:) s1(:)],-5:1:5);
    
    % Decide which is the best based on the number of near-0 intensity
    % pixels
    if hst(6,1) > hst(6,2)
        dx2.mva.(choice).img(:,:,n) = dx2.mva.(choice).img(:,:,n) * -1;
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx1,dx2] = kMeansImageCheck(dx1,dx2,choice)
% Similar idea to the PCA function above, but needs to look at the clusters
% individually

% Short cuts to the images
im1 = dx1.mva.(choice).img;
im2 = dx2.mva.(choice).img;

% How many images to do this for...
numI = size(im1,3);
    
% Loop through each image
for n = 1:numI
    
    % How many levels (clusters) are there in this image?
    numC = max(max(im1(:,:,n)));
    
    crs  = zeros(numC,numC);
    % So for each cluster
    for r = 1:numC
        
        % Extract the main image
        x = double(im1(:,:,n) == r);
        y = im2(:,:,n);
        
        
        for p = 1:numC
            
            yy = double(y == p);
            crs(p,r) = corr(x(:),yy(:));
        end
    end
    
    % So now with each cluster we find the maximum value
    maxVal = max(crs,[],1);
    [~,idx] = sort(maxVal,'descend');
        
    % What values should be changed?
    %crs2 = crs;
    map = zeros(numC,2);
    for r = 1:numC
        
        % Find the row of the largest value in column idx(r)
        [~,b] = max(crs(:,idx(r)));
        
        % Apply to the map
        map(r,:) = [b idx(r)];
        
        crs(b,:) = NaN;
                
    end
    
    % Now apply the map to the image...
    new2 = im2(:,:,n);
    for r = 1:numC
        fx = im2(:,:,n) == map(r,1);
        new2(fx) = map(r,2);
    end
    
    im2(:,:,n) = new2;

end

dx2.mva.(choice).img = im2;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resultsDisplayOptions(dpn,seg,choice)
% Presentation of the results is enacted in here.

% Just a little flourish
switch choice
    case 'PCA'
        txtLab = 'PC';
    case 'MMC'
        txtLab = 'LV';
    case 'kMeans'
        %txtLab = 'k =';
end

% Here we want to update the boxes in the panel with the number of
% components (etc) that were calculated, and a 'NA' box at the top to
% ignore the channel/colour
sz = size(dpn.d1.mva.(choice).img,3);
if strcmp(choice,'kMeans')
    labs = dpn.d1.mva.(choice).label;
    vis = {'on','off','off'};
    st = 1;
else
    
    labs = cell(sz+1,1);
    labs{1,1} = 'N/A';
    for n = 1:sz
        labs{n+1,1} = [txtLab ' ' int2str(n)];
    end
    vis = {'on','on','on'};
    st = 2;
end

% Set the values of the RGB chanel boxes
set(seg.res,'String',labs);
set(seg.res(1),'Value',min([st   sz]),'Visible',vis{1});
set(seg.res(2),'Value',min([st+1 sz]),'Visible',vis{2});
set(seg.res(3),'Value',min([st+2 sz]),'Visible',vis{3});

% Now look to refresh the images on display...
if strcmp(dpn.mode,'dual')
    dpnIonImage([],[],dpn.fig.ax.ms1,dpn.d1.mva.(choice).img(:,:,min([1 sz])));
    dpnIonImage([],[],dpn.fig.ax.ms2,dpn.d2.mva.(choice).img(:,:,min([2 sz])));
    dpnIonImage([],[],dpn.fig.ax.fu, dpn.fuse.mva.(choice).img(:,:,min([3 sz])));

%     cols = [241 182 218; 253 174 97; 26 150 65]/256;
% 
%     im1 = dpn.d1.mva.(choice).img(:,:,1:3);
%     im2 = dpn.d2.mva.(choice).img(:,:,1:3);
%     im3 = dpn.fuse.mva.(choice).img(:,:,1:3);
% 
%     im1 = predProbs2Images(im1,cols,0.95);
%     im2 = predProbs2Images(im2,cols,0.95);    
%     im3 = predProbs2Images(im3,cols,0.95);
% 
%     dpnIonImage([],[],dpn.fig.ax.ms1,im1);
%     dpnIonImage([],[],dpn.fig.ax.ms2,im2);
%     dpnIonImage([],[],dpn.fig.ax.fu, im3);

end

% Set the callback function for the three RGB boxes
set(seg.res,'Callback',{@dpnChangeSegImage2,dpn.fig,seg,choice});
dpnChangeSegImage2([],[],dpn.fig,seg,choice);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
