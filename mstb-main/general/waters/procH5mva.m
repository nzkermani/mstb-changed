function [ mva ] = procH5mva(gc,method)
%procH5pca - perform mva on the data. May need to read in from files if not
%performed directly

% Define the m/z range
mzRange = [100 1000];

% How many files?
numF = size(gc.fn,1);

% How to normalise the data?
[gc,mzMask] = doNormETC(gc,mzRange,'mean div');
[gc,~] = doNormETC(gc,mzRange,'log');

% Which analysis are we performing?
tt = tic;
switch method    
    
    case 'kmeans'
        % Standard kmeans clustering
        [mva] = km(gc);
        
    case 'pca-in'
        % In memory PCA
        [mva] = inPCA(gc);
        
    case 'pca-out'
        % Out of memory PCA, using Alan Race's memeffpca
        
    case 'nmf'
        % Need to find a good implementation of this
        [mva] = inNMF(gc);
        
    case 'tsne'
        [mva] = inTSNE(gc);
        
    otherwise
        disp('There is no otherwise');
end

mva.time = toc(tt);
mva.method = method;
mva.mzRange = mzRange;

try
    mva.mzMask = mzMask;
catch
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gc,mzMask] = doNormETC(gc,mzRange,method)

% Define m/z mask
mzMask = gc.mz >= min(mzRange) & gc.mz <= max(mzRange);

% Determine indices of individual files
[unq,~,ind] = unique(gc.idx);

switch method
    
    case {'mean div','median div'}
        
        % Define function for average spectrum
        if strcmp(method,'mean div')
            avgFunc = @nanmean;
        else
            avgFunc = @nanmedian;
        end
                    
        % Determine tissue/background pixels
        tobg = gc.tobg == 1;
        
        % Somewhere to save the average spectrum
        allAvg = zeros(numel(unq),size(gc.sp,2));

        for n = 1:numel(unq)
            
            % Determine indices for this file
            tmp = ind == n;
            
            % Determine tissue-only mean spectrum
            mnsp = avgFunc(gc.sp(tmp & tobg,:),1);
            mnsp(mnsp == 0) = 1;
                        
            % Divide all pixels by mean spectrum
            gc.sp(tmp,:) = bsxfun(@rdivide,gc.sp(tmp,:),mnsp);
            allAvg(n,:) = avgFunc(gc.sp(tmp,:),1);
            
        end
        
        % Just trim the variables here
        gc.sp = gc.sp(:,mzMask);
        allAvg = allAvg(:,mzMask);
        
        % Now determine PQN scaling factors for each of the sections based
        % on the average spectrum from each
        allFC = bsxfun(@rdivide,allAvg,avgFunc(allAvg,1));
        scFac = nanmedian(allFC,2);
        
        % Now scale each section with its scaling factor
        for n = 1:numel(unq)
            tmp = ind == n;
            gc.sp(tmp,:) = gc.sp(tmp,:) / scFac(n);
        end
        
        
%         % Now determine average tissue spectrum over all pixels... then we
%         % will do another normalisation
%         allAvg = avgFunc(gc.sp(tobg,:),1);
%         
%         % Divide to determine fold changes
%         allFC = bsxfun(@rdivide,gc.sp,allAvg);
%         
%         % Determine median fold change
%         medFC = nanmedian(allFC,2);
%         medFC(medFC == 0) = 1;
%         
%         % Finally divide
%         gc.sp = bsxfun(@rdivide,gc.sp,medFC);
        
        
    case 'tic'
        
        tt = nansum(gc.sp(:,mzMask),2);        
        gc.sp = bsxfun(@rdivide,gc.sp(:,mzMask),tt) * 1000;
        
    case 'log'
        
        os = nanmedian(gc.sp(gc.sp > 0));        
        gc.sp = log(gc.sp + os);
        
    otherwise
        disp('No otherwise');
        
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mva] = inPCA(gc)
% In memory PCA

% Apply the mask to the background pixels
mask = gc.tobg == 1;

[mva.ll,mva.ss,mva.ee] = pca(gc.sp(mask,:),...
    'NumComponents',10,...
    'Economy',true);

% Need to reconstruct after PCA using only TO (not BG) pixels...
tmp = zeros(size(gc.idx,1),size(mva.ss,2));
tmp(mask,:) = mva.ss;
mva.ss = tmp;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mva] = inNMF(gc)
% In memory NMF

% Apply the mask to the background pixels
mask = gc.tobg == 1;

% Initial random NMF
t1 = tic;
opts = statset('Display','final',...
    'UseParallel',true,...
    'MaxIter',10);
[initSS,initLL] = nnmf(gc.sp(mask,:),5,...
    'Options',opts,...
    'Replicates',20,...
    'Algorithm','mult');
t1 = toc(t1);

% Do a more accurate version
t2 = tic;
opts = statset('Display','final',...
    'UseParallel',true,...
    'MaxIter',20);
[mva.ss,mva.ll,mva.resid] = nnmf(gc.sp(mask,:),5,...
    'w0',initSS,...
    'h0',initLL,...
    'Options',opts,...
    'Replicates',10,...
    'Algorithm','als');
t2 = toc(t2);

[t1 t2]

% Need to reconstruct after PCA using only TO (not BG) pixels...
tmp = zeros(size(gc.idx,1),size(mva.ss,2));
tmp(mask,:) = mva.ss;
mva.ss = tmp;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mva] = km(gc)
% In memory PCA

% Only run on tissue pixels
mask = gc.tobg == 1;

% Run kmeans
idx = kmeans(gc.sp(mask,:),4);

% Reconstruct...
tmp = zeros(size(gc.idx,1),1);
tmp(mask,:) = idx;
mva.ss = tmp;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mva] = inTSNE(gc)

mva = [];
return

% Only run on tissue pixels
mask = gc.tobg == 1;

% TSNE here...
    
% Parameters... need to work out what these actually do
perplexity = 30;
layers = [500 500 2000 2];
    
% Determine test and train indices...
xx = [];

% Train the parametric t-SNE network
[network, err] = train_par_tsne(train_X, train_labels, test_X, test_labels, layers, 'CD1');


% Reconstruct...
tmp = zeros(size(gc.idx,1),1);
tmp(mask,:) = idx;
mva.ss = tmp;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

