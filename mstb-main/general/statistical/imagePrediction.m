function [train,test,mmc] = imagePrediction( train,test )
%imagePrediction - use one MS image to predict other ones.
%
% INPUTs
% Two structures, train and test, which contain the following fields:
% name      - file name
% mz        - [1 x n] vector of n m/z values
% sp        - [p x q x n] array of unlogged spectral intensities
% *anno     - [p x q x c] array of annotations for c classes
% *groups   - {1 x c} cell array of c class labels
%
% The latter two fields are required only for the train variable.
% The test structure can contain multiple levels, i.e. test(n).sp
%
% OUTPUTs
% Some images probably
%
% DATA
% Exemplar data can be found at:
% /Users/jmckenzi/Google Drive/Jocelyn/RatBrain.mat
%
% USAGE
% imagePrediction(train,test);
%
% James McKenzie, 2016

% Define parameters for the entire function, which could be made better
% with a more refined function
opts.ppmTol = 10;
opts.numComps = 3;
opts.mzRange = [0 Inf];
opts.doLog   = true;

% Let's set the background pixels to NaN values and totally ignore them
train.tobg = pixTOBG(nansum(train.sp,3),[],true);
for n = 1:size(test,2)
    test(n).tobg  = pixTOBG(nansum(test(n).sp,3),[],true);
end

% Actually first we should consider reshaping the matrices
[train] = reshapeMatrices(train);
[test]  = reshapeMatrices(test);

% First we need to do the peak matching between the m/z vectors
[test,train] = mzAlignment(train,test,opts.ppmTol,opts.mzRange);

% Run the MMC model
[mmc.B,mmc.W,mmc.trainMean] = mmcModel(train,opts);
mmc.opts = opts;
% Now predict a model to the things in the other files...

return

%return

% Perform MMC on this newly aligned data
[train,test] = prepMMC(train,test,opts);

% Generate the images and hey presto!
[train] = reSquareMatrices(train);
[test]  = reSquareMatrices(test);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x] = reshapeMatrices(x)
% Reshape from images to matrices

for n = 1:size(x,2)
    x(n).sz = size(x(n).sp);

    x(n).sp = reshape(x(n).sp,[x(n).sz(1)*x(n).sz(2) x(n).sz(3)]);

    if ~isempty(x(n).anno)
        x(n).anno = reshape(x(n).anno,[x(n).sz(1)*x(n).sz(2) size(x(n).anno,3)]);
    end
    
    if isfield(x(n),'tobg')
        x(n).tobg = reshape(x(n).tobg,[x(n).sz(1)*x(n).sz(2) 1]);
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diagnosticSpectra(train,test)

cols = jet(size(test,2) + 1);

figure; hold on;

% Find non-bg pixels and calculate the sum
sp0 = nanmean(train.sp(train.tobg,:),1);
plot(train.mz,sp0,'-o','Color',cols(1,:));

for n = 1:size(test,2)
    plot(test(n).mz,nanmean(test(n).sp(test(n).tobg,:),1),'-o',...
        'Color',cols(n+1,:));
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diagnosticSpectraAL(train,test)

cols = jet(size(test,2) + 1);

figure; hold on;

% Find non-bg pixels and calculate the sum
sp0 = nanmean(train.sp(train.tobg,:),1);
plot(train.mz,sp0,'-o','Color',cols(1,:));

for n = 1:size(test,2)
    plot(test(n).mzAl,nanmean(test(n).spAl(test(n).tobg,:),1),'-o',...
        'Color',cols(n+1,:));
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [test,train] = mzAlignment(train,test,ppmTol,mzRange)
% There are two potential schools of thought regarding m/z alignment.  We
% can either go for the pre-defined cmz as the training set m/z vector, or
% find the consensus between the samples and align everything in a common
% fashion.  I'll implement the first method initially, but hope to have
% both done eventually

% Define the peak matching options
opts.handBand = eval(['@(z) (' num2str(ppmTol) '* z / 1e6)']);            
opts.handGap  = eval(['@(z) (' num2str(ppmTol) '* z(:,1) / 1e6)']);
opts.handDist = @(R,S) abs(sum((R-S),2));
opts.boolSC = false; 
opts.boolSA = false; 
opts.boolSN = false;

mask = train.mz >= min(mzRange) & train.mz <= max(mzRange);

% This is the required format for the cmz vector and intentity
cmz = [train.mz(mask)' max(full(train.sp(:,mask)),[],1)'];

% The peak alignment procedure essentially deletes peaks from others, but
% essentially removes their intensities from the files, and these were
% important contributison when they were originally normalised.  So need to
% re-scale the files a little bit
numF = size(test,2);
allMean = zeros(numF+1,size(cmz,1));
allMean(1,:) = nanmean(train.sp(:,mask),1);

figure; hold on;
cols = jet(numF);
% Loop through each file in turn and align to the cmz
for n = 1:numF
    
    stem(test(n).mz,nanmean(test(n).sp,1),'Color',cols(n,:));
   
    % This is the mz and intensity vector of sample n
    tmpMZ = [test(n).mz' max(full(test(n).sp),[],1)'];
    
    [j,k] = samplealign2(...
        cmz,...
        tmpMZ,...
        'Band',opts.handBand,...
        'Gap',opts.handGap,...
        'Distance',opts.handDist,...
        'Quantile',[],...
        'SHOWCONSTRAINTS',opts.boolSC,...
        'SHOWNETWORK',opts.boolSN,...
        'SHOWALIGNMENT',opts.boolSA);

    % Make the matrices to finish...
    newSP = zeros(size(test(n).sp,1),size(cmz,1));
    newSP(:,j) = test(n).sp(:,k);
    
    % Save the new matrix and new mz vector
    test(n).mzAl = train.mz(mask);
    test(n).spAl = newSP;
    test(n).szAl = test(n).sz;
    test(n).szAl(3) = sum(mask);%size(train.sp,2);
    
    % Determine the mean spectrum here...
    allMean(n+1,:) = nanmean(newSP,1);
       
    
end
stem(train.mz,nanmean(train.sp,1),'Color','k');

[~,scaleFac] = jsmNormalise(allMean(2:end,:),'pqn-mean',0,0,[],allMean(1,:));

% Normalise to this scaling factor in order to get everything on the same
% approximate scale
for n = 1:numF
    test(n).spAl = test(n).spAl ./ scaleFac(n);
    test(n).newlogOS = test(n).logOS ./ scaleFac(n);
end

% Update sizes as dependent on the m/z range used...
train.sp = train.sp(:,mask);
train.mz = train.mz(mask);
train.sz(3) = sum(mask);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B,W,meanTrain] = mmcModel(train,opts)
% Prepare the data for MMC and then run it using the annotated pixels

% Generate the matrix of [0 1 2 n] which denotes the class of the pixel,
% where 0 is for unannotated pixels which are to be predicted
trainID = bsxfun(@times,train.anno,1:size(train.anno,2));
trainID = sum(trainID,2);

% Determine indices of the training set
idxTrain = trainID ~= 0;

% Determine the actual training set
spTrain = train.sp(idxTrain,:);

% Log transformation?
if opts.doLog
    spTrain = log(spTrain + train.logOS);
end

% What is the mean spectrum?
meanTrain = nanmean(spTrain,1);

% Subtract the mean spectrum from the training set
spTrain = bsxfun(@minus,spTrain,meanTrain);

% Quick and dirty PCA
qdPCA(spTrain,trainID(idxTrain));

% Run the OAA-MMC function
[B,W] = oaaMMC(spTrain,trainID(idxTrain),opts);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [train,test] = prepMMC(train,test,opts)
% Prepare the data for MMC and then run it using the annotated pixels

% Generate the matrix of [0 1 2 n] which denotes the class of the pixel,
% where 0 is for unannotated pixels which are to be predicted
trainID = bsxfun(@times,train.anno,1:size(train.anno,2));
trainID = sum(trainID,2);

% Loop through the test set, generating sample IDs for each
% Also determine a more suitable log offset for each file type
numF = size(test,2);
fInfo = struct('sampleID',[],'logOS',[]);
for n = 1:numF
    fInfo(n).sampleID = ones(test(n).sz(1)*test(n).sz(2),1) * n;
    
    % Replicate the newlogOS 
    fInfo(n).logOS = repmat(test(n).newlogOS,[test(n).sz(1)*test(n).sz(2),1]);
end
sampleID = vertcat(fInfo.sampleID);

% Concatenate all of the other pixels in test into a mega matrix
allTest = vertcat(test.spAl);
testID = zeros(size(allTest,1),1);

% Prepare the mega matrix
allSpec = [train.sp; allTest];
allHist = [trainID; testID];
allSamp = [zeros(size(train.sp,1),1); sampleID];
%allLogs = [repmat(train.logOS,[size(train.sp,1) 1]); newlogOS];

% Determine indices of the training set
idxTrain = allHist ~= 0;%& nansum(allSpec,2) > 0;
%meanTrain = nanmean(allSpec(idxTrain,:),1);
meanTrain = nanmean(allSpec(:,:),1);

% Subtract the mean spectrum from the training set
spTrain = bsxfun(@minus,allSpec(idxTrain,:),meanTrain);

% Quick and dirty PCA
qdPCA(train,trainID);

% Run the OAA-MMC function
[B,W] = oaaMMC(spTrain,idxTrain,allHist,opts);



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doMCM
% Calculate probabilities for pixels
[T,Pr,~] = calcScoresSpectra(allSpec,nanmean(meanTrain,1),W,B);





% So now that we have the probabilities for each of the sections in this
% mega matrix, we need to split it up again according to sampleID
for n = 0:1:numF
    
    % Indices...
    fx = allSamp == n;
    
    if n == 0
        train.prob = Pr(fx,:);
        train.score = T(fx,:);
    else
        test(n).prob = Pr(fx,:);
        test(n).score = T(fx,:);
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function qdPCA(train,trainID)
% QUICK and dirty PCA of just the training set's annotated pixels

[~,ss,~] = princomp(train,'econ');

figure; hold on;

scatter3(ss(:,1),ss(:,2),ss(:,3),80,trainID,'o','filled');


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B,W] = oaaMMC(spTrain,groupID,opts)

% Now we'll begin the OAA approach
[unq,~,ind] = unique(groupID);
numG = numel(unq);

% Somewhere for BETA coefficients
B = zeros(2,numG);       
    
% Matrix for weights
W = zeros(size(spTrain,2),numG);

% Loop
for n = 1:numG
    
    % This is the group that is going up against the others
    tmpIdx = ind == n;
    
    % Run the function...
    [~,T,~,Wi] = recursiveMmcLda(spTrain,tmpIdx,opts.numComps);

    % Single function for the regression coefficients
    [B,W,~] = getRegCoef(n,T,tmpIdx,W,Wi,spTrain,B);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x] = reSquareMatrices(x)
% Convert the linear data matrices back into an image format

numF = size(x,2);

for n = 1:numF
    
    sz = x(n).sz;
    
    % These are common to all of the files
    x(n).sp = reshape(x(n).sp,sz);    
    x(n).prob = reshape(x(n).prob,[sz(1) sz(2) size(x(n).prob,2)]);
    x(n).score= reshape(x(n).score,[sz(1) sz(2) size(x(n).score,2)]);
    
    % These are only in the test sets
    if isempty(x(n).anno)        
        sz2 = x(n).szAl;        
        x(n).spAl = reshape(x(n).spAl,sz2);
    end    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


