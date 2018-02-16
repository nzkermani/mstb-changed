function [ data ] = predictionMMC(data,opts)
%predictionMMC - does the donkey work

numF = size(data,2);

% Create a structure to contain information regarding the file
info = struct('file',[],'bg',[],'anno',[]);
for n = 1:numF    
    sz = size(data(n).spAl,1);
    info(n).file = ones(sz,1) * n;
    
    if size(data(n).tobg,2 > 1)
        info(n).bg = reshape(data(n).tobg,[prod(data(n).sz(1:2)) 1]);
    else
        info(n).bg = data(n).tobg;
    end
    
    if ~isempty(data(n).anno)
        tmp = bsxfun(@times,data(n).anno,1:size(data(n).anno,2));
        tmp = sum(tmp,2);        
        info(n).anno = tmp;
    else
        info(n).anno = zeros(size(data(n).spAl,1),1);
    end
end

% Concatenate all data and information
allD = vertcat(data.spAl);
file = vertcat(info.file);
anno = vertcat(info.anno);
tobg = vertcat(info.bg);


% Set alll bg pixels to 0
allD(~tobg,:) = NaN;
anno(~tobg,:) = 0;

% Remove NaN values from the data
allD(isnan(allD)) = 0;

% Log the data?
if opts.doLog
    logOS = nanmedian(allD(allD > 0));
    allD = log(allD + logOS);
end

% Indices of the training set
idxTrain = anno ~= 0;

% What is the mean spectrum?
meanTrain = nanmean(allD(idxTrain,:),1);

% Run the OAA-MMC function
[B,W] = oaaMMC(allD(idxTrain,:),anno(idxTrain),opts.numComps);

% Calculate probabilities for pixels
[T,Pr,~] = calcScoresSpectra(allD,meanTrain,W,B);

% Separate according the files
for n = 1:numF
    
    % Indices...
    fx = file == n;
    
    % Size to which to reshape?
    sz = data(n).sz;
    sz(3) = size(Pr,2);
    
    data(n).prob  = reshape(Pr(fx,:),sz);
    data(n).score = reshape(T(fx,:),sz);
    
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


