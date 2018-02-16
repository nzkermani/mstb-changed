function h5Extract(fold,name,ions,opts)
%h5Extract - extract multiple ion images from an H5 file.  Display the
%optical image in conjunction with these extracted images.
%
% INPUTs
% fold - h5 location 
% name - h5 file name
% ions - vector of ions to be extracted
% opts - structure containing various options
%           .ppm - ppm tolerance

% These are the default options
if isempty(opts)
    opts.ppm = 20;
end
if ~isfield(opts,'ppm')
    opts.ppm = 20;
end


% This gets the actual images
[img,opt,anno,grps,probs,extrMZ] = getMZvectors(fold,name,ions,opts);

% Determine the values for each of the annotated regions and the pixel
% predicted regions too...
[res,heads] = annotatedIntensity(img,anno,grps,probs);

% Now we need to work on a way to display all the images...
plotImages(img,opt,anno,grps,ions,opts.ppm,probs,extrMZ,res,heads)

% Decide if we are going to save it...
if isfield(opts,'save')
    
    % Either provided a full name or a folder
    if exist(opts.save,'dir')
        saveName = [opts.save name(1:end-3)];
    else
        saveName = opts.save;
    end

    % Save the file
    graphFormat(saveName,'png');    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [img,opt,anno,grps,probs,extrMZ] = getMZvectors(fP,fN,ions,opts)
% Open the h5 files and get the mz vectors of each of the pixels
    
% File name
tmp = [fP fN];
disp(tmp);

% Optical image please
opt = h5read(tmp,'/opimage_aligned');

% Read in the m/z vector
mz = h5read(tmp,'/mz');

% Read in the full data matrix - unlog as necessary
data = h5read(tmp,'/X');
if max(data(:)) < 500
    os = min(data(:));    
    data = exp(data) - exp(os);    
end

% Perform TIC normalisation over the sum of each pixel
%ticSum = nansum(data,3);
%data = bsxfun(@rdivide,data,ticSum) * 1e6;

% Determine the size of the image
sz = size(data);

% Somewhere to store the data
numI = numel(ions);
if numI == 1
    %img = zeros([sz(1) sz(2)]);
else
    img = zeros([sz(1) sz(2) numI]);
end

% Somewhere to save the mz values that we actually extract for each image
extrMZ = cell(numI,1);

% Now loop through the mz values that we are supposed to be extracting
for n = 1:numI

    % Ion tolerances
    il = ions(n) * (1e6 - opts.ppm)/1e6;
    ih = ions(n) * (1e6 + opts.ppm)/1e6;

    % Gather the correct ions
    [fx] = mz > il & mz < ih;
    
    if sum(fx) == 0
        continue;
    else
        extrMZ{n,1} = mz(fx);
    end        

    % Place into the images
    if numI == 1
        img = sum(data(:,:,fx),3);
    else
        img(:,:,n) = sum(data(:,:,fx),3);
    end            
end
    
% This is the total ion image which is useful for the determination
% of the tissue object pixels
fx = mz > 885.45 & mz < 885.65;
tot = sum(data(:,:,fx),3);

% Now here we do the annotations - which gets complicated when we need to
% compare the files
try
    anno = h5read(tmp,'/groupPixels');
catch
    anno = zeros(sz(1),sz(2));
    grps = {'None'};
    return
end

% Now the actual annotation names
numA = size(anno,3);
grps = cell(numA,1);
for n = 1:numA
    grps{n,1} = h5readatt(tmp,'/tissue_id',int2str(n));
end

% What about the pixel probabilities based on MMC predictions?
try
    probs = h5read(tmp,'/pixelProbs');
catch
    disp('probs failed');
    probs = zeros(sz(1),sz(2),numA);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pixDiag,intVals] = detPixGroups(img)
% For each of the ion images, determine the connectivity of the pixels and
% delete pixels which aren't connected to minPix neighbours
% Perform the connectivity function to ensure that we have at least two
% neighbouring pixels...

% Loop through each image
numI = size(img,3);

pixDiag = zeros(numI,2);
intVals = zeros(numI,3);

% Create a filter of [3 3] size for smoothing out the image
filt = fspecial('average',[3 3]);

for n = 1:numI
    
    % Just look at non-zero pixels
    tmp = img(:,:,n) > 0;
    
    % Now apply the filtering
    f2 = imfilter(tmp,filt);
    
    % Determine the intensities of the TO pixels
    tmp2 = img(:,:,n);
    tmp2 = tmp2(f2);
    intVals(n,:) = prctile(tmp2(~isnan(tmp2)),[25 50 75]);
    
    % Now do the pixel connectivity
    conn = bwconncomp(f2);
    
    % These are the sizes of the connected pixels
    szes = cellfun(@size,conn.PixelIdxList','UniformOutput',false);
    
    % See what we find
    if ~isempty(szes)
        
        % These are the sizes of the islands
        szes = max(vertcat(szes{:}),[],2);
        
        pixDiag(n,1) = max(szes(:));
        pixDiag(n,2) = sum(szes);                
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res,heads] = annotatedIntensity(img,anno,grps,prob)
% Determine the intensities of the annotated regions

thresh = 0.99;

sz = size(img);
numA = numel(grps);
numI = sz(3);

% Start by reshaping everything
img  = reshape(img ,[sz(1)*sz(2) numI]);
anno = reshape(anno,[sz(1)*sz(2) numA]);
prob = reshape(prob,[sz(1)*sz(2) numA]);

% Convert the pixel probs so that there are no ambiguous pixels included
mask = prob >= thresh;
tmp  = sum(mask,2) > 1;
mask(tmp,:) = false;

% Store all of the data
res = NaN(numI,10,numA);

for n = 1:numA
    
    % Annotated-only pixels
    %grps{n}
    fx = anno(:,n) == 1;
    
    vals = img(fx,:);    
    med1 = prctile(vals,[25 50 75],1)';
    mnn1 = mean(vals,1);
    std1 = std(vals,[],1);
    
    
    % What about the predicted pixels
    vals2 = img(mask(:,n),:);
    med2 = prctile(vals2,[25 50 75],1)';
    mnn2  = mean(vals2,1);
    std2  = std(vals2,[],1);
    
    % Save the data...
    res(:,:,n) = [med1 mnn1' std1' med2 mnn2' std2'];
    
end
    
% Headers so that we remember what the numbers are later on;
heads = {'Anno25','Anno50','Anno75','AnnoMean','AnnoSTD',...
    'Prob25','Prob50','Prob75','ProbMean','ProbSTD'};
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = plotImages(img,opt,anno,grps,ions,ppmTol,prob,extrMZ,res,heads)
% Mega plotting function...

fig.fig = figure('Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5],...
    'Name','H5 EXTRACT',...
    'Number','off');

fig.pan1 = uipanel('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0 0 0.25 1],...
    'BackgroundColor','black');

fig.pan2 = uipanel('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.25 0 0.75 1],...
    'BackgroundColor','black');

% Two axes for the optical images in panel 1
fig.op1 = axes('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0 0.5 1 0.5]);
imagesc(opt);
axis square

% Now the annotation regions - needs to be prettified
fig.op2 = axes('Parent',fig.pan1,...
    'Units','normalized',...
    'Position',[0 0 1 0.5]);
anno2 = bsxfun(@times,anno,reshape(1:size(anno,3),[1 1 size(anno,3)]));
anno2 = nansum(anno2,3);
imagesc(anno2);
axis square

% Determine arrangement for image grid
numI = size(img,3);
optRC = sqrt(numI);
if isinteger(optRC)
    arr = [optRC optRC];
else
    tmp2 = ceil(optRC);
    optRC = ceil(numI / tmp2);
    arr = [optRC tmp2];
end


% Draw each ion
spax = zeros(numI,1);
for n = 1:numI
    
    % Subplot here
    spax(n,1) = subplot(arr(1),arr(2),n,...
        'Parent',fig.pan2);
        
    % Determine maximum value
    maxI = max(max(img(:,:,n)));

    % This actually plots the ions...
    if maxI > 0
        imagesc(img(:,:,n));
    else
        imagesc(NaN(size(img,1),size(img,2),3));
    end
    axis square
    
%     % Global colour limit application
%     if cLim
%         caxis(colLims);
%     end
    
    % Need to have a title for each of the images in order to know what it
    % is showing us.
    ppmDiff = ppmTol * ions(n) / 1e6;
    titTxt = ['m/z = ' sprintf('%0.3f',ions(n)) ' ± ' ...
        sprintf('%0.3f',ppmDiff)];
    title(titTxt,...
        'FontSize',16,...
        'Color',[0.5 0.5 0.5],...
        'FontWeight','bold');
    
    text(2,size(img,1)-1,sprintf('%0.0f',maxI),...
        'Color',[0.5 0.5 0.5]);
    
end
colormap(hot);

% Text string for histo classes
allCols = hot(numel(grps));
for n = 1:numel(grps)

    % Determine defaults...
    if n == 1 
        colLab = ['\fontsize{16} {\color[rgb]{' sprintf('%f %f %f',allCols(n,:)) '}' grps{n} ' }'];
    else
        colLab = [colLab ' | {\color[rgb]{' sprintf('%f %f %f',allCols(n,:)) '}' grps{n} ' }'];
    end
end

fig.axTxt = axes('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0 0 0.001 0.03],...
    'Color','black');
text(1,1,colLab,...
    'FontWeight','bold');

% Axes formatting
set([fig.op1 fig.op2 spax'],...
    'XTick',[],...
    'YTick',[],...
    'XColor',[0.2 0.2 0.2],...
    'YColor',[0.2 0.2 0.2]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



