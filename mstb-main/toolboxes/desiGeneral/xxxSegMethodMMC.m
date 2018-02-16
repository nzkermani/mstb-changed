function [dx] = xxxSegMethodMMC(dx,wb,anno)
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
    if ~isempty(wb)
        waitbar(n/numD,wb,['MMC: ' int2str(n) '/' int2str(numD)]);
    end
    
    % Information about this classification
    hh = hist(mmcNumbs,unique(mmcNumbs))
    sum(hh)
    100 * hh / sum(hh)

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

