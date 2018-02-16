function [pq,hi] = statsUnivariateMethod(mz,sp,histID,method)
% Calculate univariate statistics on the data

% Somewhere to store the results
numV = size(sp,2);
pq = zeros(numV,2);

wb = waitbar(0,'Univariate calculations!');

warning off all

% Run the various methods
switch method    
    case {'ANOVA'}
        for n = 1:numV
            pq(n,1) = anova1(sp(:,n),histID,'off');
            waitbar(n/numV,wb);
        end
        
    case {'Kruskal-Wallis'}
        for n = 1:numV
            pq(n,1) = kruskalwallis(sp(:,n),histID,'off');
            waitbar(n/numV,wb);
        end
        
    case 'Output'
        pqa = zeros(numV,2);
        pqk = zeros(numV,2);
        auc = zeros(numV,1);
        [~,~,unqInd] = unique(histID);
        for n = 1:numV
            pqa(n,1) = anova1(sp(:,n),histID,'off');
            pqk(n,1) = kruskalwallis(sp(:,n),histID,'off');
            [~,~,~,auc(n,1)] = perfcurve(unqInd,sp(:,n),2);
            if auc(n,1) < 0.5
                auc(n,1) = 1 - auc(n,1);
            end
            waitbar(n/numV,wb);
        end
        
        % Now determine the q values
        [pqa(:,2),~] = getBHYqVls(pqa(:,1)',0.05);
        [pqk(:,2),~] = getBHYqVls(pqk(:,1)',0.05);
        
        % Dump these variables to the workspace
        dump = [pqa pqk auc];
        
        assignin('base','dump',dump);
        
        pq = pqa;
        
    otherwise
        error('there is no otherwise');
end

warning on all

% Perform FDR on the p values
[pq(:,2),~] = getBHYqVls(pq(:,1)',0.05);

% Whilst there are possibly multiple groups in this analysis, let's try to
% calculate a fold change of some description
[unq,~,ind] = unique(histID);
numG = numel(unq);
fc = zeros(numG,numV);
for n = 1:numG
    
    % Mean spectrum of reference group
    fx = ind == n;
    ms = nanmean(sp(fx,:),1);
    
    % Store the FCs here
    tmpFC = NaN(numG,numV);
    
    for r = 1:numG        
        if n ~= r
            
            % Mean spectrum of this group
            fy = ind == r;
            rf = nanmean(sp(fy,:),1);
            
            % Calcualte the log2 fold changes
            tmpFC(r,:) = log2(ms ./ rf);
        end
    end
    
    % Remove NaN / Inf
    tmpFC(isnan(tmpFC)) = 0;
    tmpFC(isinf(tmpFC)) = 0;
    
    % Find largest/smallest
    lo = nanmin(tmpFC,[],1);
    hi = nanmax(tmpFC,[],1);    
    fx = abs(lo) > hi;    
    hi(fx) = lo(fx);
    
    % These are the biggest FCs cf this group
    fc(n,:) = hi;
    
    
end
       
lo = nanmin(tmpFC,[],1);
hi = nanmax(tmpFC,[],1);
fx = abs(lo) > hi;
hi(fx) = lo(fx);
           
% Delete the waitbar
delete(wb);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
