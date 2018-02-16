function [scfac] = procH5normFactor(file)
%procH5normFactor - read in the raw data sequentially in order to determine
%suitable normalisation / transformation scaling values...

% File name
if nargin == 0
    %file = '/Users/jmckenzi/Documents/Box Sync/H5 Proc/H5-Test-170510-135411.h5';
    file = '/Users/jmckenzi/Documents/Box Sync/H5 Proc/PP-FullRange-10-Proc-170519-121108.h5';
end

refSpec = 'mean';

% Work out how many datasets are included within...
info = h5info(file);
numF = size(info.Groups,1);

% Read in the mz vector
fn = cell(numF,1);
%mz = h5read(file,'/mz');
mz = h5read(file,['/file1/rawMZ']);

% Create an m/z mask to remove anything to do with raffinose, which tends
% to distort everything quite markedly
mask = ~(mz > 500 & mz < 600);

% Create a matrix to store, e.g. mean spectra
avgSpec = zeros(numF,numel(mz));

% Read in each part...
sp = cell(numF,3);
for n = 1:numF
    
    % Read the name of the raw file
    fn{n} = h5readatt(file,['/file' int2str(n)],'fileName');
    
    % Read in the spectral matrix
    tmp = h5read(file,['/file' int2str(n) '/img']);
    sz = size(tmp);
        
    % Reshape it to normal [m x n] matrix
    tmp = reshape(tmp,[sz(1)*sz(2) sz(3)]);
        
    % Determine the tissue and background pixels. Use 2-means clustering
    % and pick the cluster with the highest intensity
    %[tobg] = determineBackground(tmp,mask,sz);

    % Do something depending on the reference required
    switch refSpec
        case 'mean'
            avgSpec(n,:) = nanmean(tmp(:,:),1);
        case 'median'
            avgSpec(n,:) = nanmedian(tmp(:,:),1);
    end
    
end

% Determine the PQN scaling factor
fc = bsxfun(@rdivide,avgSpec,nanmean(avgSpec,1));
scfac = nanmedian(fc,2);

% scaled = bsxfun(@rdivide,avgSpec,scfac);
% 
% [ll1,ss1,ee1] = pca(avgSpec,'NumComponents',2);
% ee1 = 100 * ee1 / sum(ee1);
% 
% [ll2,ss2,ee2] = pca(scaled,'NumComponents',2);
% ee2 = 100 * ee2 / sum(ee2);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tobg] = determineBackground(img,mask,sz)
% Determine the background pixels using 2 means clustering

[a] = dpnTOBG(nansum(log(img(:,mask)+1),2),[],[]);

idx = kmeans(log(img(:,mask)+1),2);

% Summed spectra from each cluster
k1 = sum(img(idx == 1,mask),1);
k2 = sum(img(idx == 2,mask),1);

% Determine fold changes
fc = log2(k1 ./ k2);

% Mean pos/neg values
mpos = abs(mean(fc(fc > 0)))
mneg = abs(mean(fc(fc < 0)))

if mpos > mneg
    tobg = idx == 1;
else
    tobg = idx == 2;
end

hh = reshape(tobg,[sz(1) sz(2)]);
aa = reshape(a,[sz(1) sz(2)]);

figure; imagesc(hh);
figure; imagesc(aa);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%