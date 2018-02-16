function predictionPCA(data,opts)
%predictionPCA - run PCA on the newly m/z aligned data in order to see if
%it needs to be normalised etc...
%
% James McKenzie

numF = size(data,2);

% Create a structure to contain information regarding the file
info = struct('file',[],'bg',[]);
for n = 1:numF    
    sz = size(data(n).spAl,1);
    info(n).file = ones(sz,1) * n;    
    info(n).bg = data(n).tobg;    
end

% Concatenate all data and information
allD = vertcat(data.spAl);
file = vertcat(info.file);
tobg = vertcat(info.bg);
inds = 1:20:size(allD,1);

% Set alll bg pixels to 0
allD(~tobg,:) = NaN;

% Log the data?
if opts.doLog
    logOS = nanmedian(allD(allD > 0));
    allD = log(allD + logOS);
end

% This is the mean spectrum
%meanSp = nanmean(allD(inds,:),1);

% Run the PCA
[ll,~,ee] = princomp(allD(inds,:),'econ');
ss = allD * ll;

% Break the PCA up into files
figure; hold on;
cols = jet(numF);
for n = 1:numF
    
    fx = file == n;
    
    tmp = ss(fx,1:3);
    
    idx = randperm(size(tmp,1),600);
    scatter(tmp(idx,1),tmp(idx,2),80,cols(n,:),'o','filled');
    
    data(n).pca = reshape(tmp,[data(n).szAl(1) data(n).szAl(2) size(tmp,2)]);
    
end

% Now plot together
figure;
ax = zeros(numF,1);
for n = 1:numF
    
    ax(n,1) = subplot(1,numF,n);
    tmp = data(n).pca(:,:,1:3);
    tmp = imScale(tmp);
    imagesc(tmp);
    
end

figure; stem(data(1).mzAl,ll(:,1));


end

