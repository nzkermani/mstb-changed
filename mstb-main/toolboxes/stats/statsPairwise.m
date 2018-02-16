function [ res ] = statsPairwise(msa)
%statsPairwise - do this for Anna's PBC samples to compare one condition to
%each other condition.

% We only want to use the nodules, so find only observations with that word
fx = ~cellfun(@isempty,strfind(lower(msa.meta.histID),'nodules'));

% Strip out only these observations
sp = msa.sp(fx,:);
histID = msa.meta.foldID(fx,:);
fileID = msa.meta.fileID(fx,:);

% Determine unique histIDs
[unq,~,ind] = unique(histID);
numH = numel(unq);

% Confusion matrix for classification results
res = NaN(numH,numH);

% Loop through
for n = 1:numH
    
    % Indices for group n
    indN = ind == n;
    
    disp([unq{n} ' ----> ' int2str(sum(indN))]);
    
    % Begin classification against subsequent groups
    for r = n+1:numH
        
        % Indices for group r
        indR = ind == r;
        
        % Combined indices
        comb = indN | indR;
        
        % Extract spectra and file/histIDs
        tmpSp = sp(comb,:);
        tmpHist = histID(comb,:);
        tmpFile = fileID(comb,:);
        
        % Perform TIC normalisation
        tmpSp = bsxfun(@rdivide,tmpSp,nansum(tmpSp,2)) * 1000;
        os = nanmedian(tmpSp(tmpSp > 0));
        
        % Leave file out CV
        [~,~,allP,~,~] = lpoMMC(tmpSp,tmpHist,tmpFile);
        
        [~,~,allP2,~,~] = lpoMMC(log(tmpSp+os),tmpHist,tmpFile);
        
        % Now determine total classification performance
        cor = allP(:,1) == allP(:,2);
        res(r,n) = sum(cor) / numel(cor);
        
        cor2 = allP2(:,1) == allP2(:,2);
        res(n,r) = sum(cor2) / numel(cor2);
        
        
        % % Plot the scores / confusion matrix
        % [a,b,c] = unique(tmpHist);
        % 
        % figure; hold on;
        % scatter(sc(c == 1,1),sc(c == 1,2),'b','o','filled');
        % scatter(sc(c == 2,1),sc(c == 2,2),'r','o','filled');
        
        
    end
    
end
 
res = res * 100;

figure; imagesc(res);
caxis([0 100]);
cb = colorbar;

        

end

