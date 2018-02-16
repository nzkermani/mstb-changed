function [ allS,cm,allP,unq,allL ] = kfoldMMC(sp,grp,k)
%kfoldMMC

verb = false;
if verb
    f1 = figure;
    ax = axes('Parent',f1);
end

% Create a waitbar
wb = waitbar(0,'kFold CV');

% Create the CV partition data matrix
cvp = cvpartition(grp,'k',k);

% Number of spectra
[numS,numV] = size(sp);
trainVec = (1:numS)';

% How many DVs to calculate?
[unq,~,allIdx] = unique(grp);
numG = numel(unq);
numDV = max([2 numG-1]);

% Check that there are enough groups
if numG == size(sp,1)
    disp('Cannot perform supervised analysis');
    delete(wb);
    return
end

% Somewhere to store the predicted scores...
allS = zeros(numS,numDV);
allP = zeros(numS,2);
allL = zeros(numV,numDV);

% Run a simple, single analysis using ALL available observations
[~,~,~,totLL] = recursiveMmcLda(sp,allIdx,numDV);

% Confusion matrix...
cm = zeros(numG,numG);

% Loop through the folds
for n = 1:k
    
    % Define the training/testing set indices
    idxTrain = cvp.training(n);
    idxTest  = cvp.test(n);
    
    % Define the groupings
    [~,~,grpIdx] = unique(grp(idxTrain));
    
    % Run the MMC on the training set
    [~,ss,~,ll] = recursiveMmcLda(sp(idxTrain,:),grpIdx,numDV);
    
    % This function (replacing the previous one) needs to map the
    % calculated scores on to the overall scores to ensure consistency.
    % This is because axes can be flipped which will ruin the average
    % scores to be presented at the end of the analysis...
    [~,vals] = alignCVdirsDR(totLL,ll);
    ll = bsxfun(@times,ll,vals);
    ss = bsxfun(@times,ss,vals);
    
%     % We should mirror the scores in the case of 2 groups in order that
%     % they are compatible.  This prevents the scores from looking stupid if
%     % the first axes gets swapped around. I'm not sure how this will be
%     % shown if there are more than 2 groups in the analysis
%     if numG == 2
%         if n == 1
%             firstMed = median(ss(grpIdx == 1,1));
%         else
%             thisMed = median(ss(grpIdx == 1,1));
%             if (firstMed < 0 && thisMed > 0) || (firstMed > 0 && thisMed < 0)
%                 ss(:,1) = ss(:,1) * -1;
%                 ll(:,1) = ll(:,1) * -1;
%             end
%         end
%     end

    % Save the loadings
    allL = allL + ll;
    
    % Determine the mean spectrum of the training set
    meanTrain = nanmean(sp(idxTrain,:),1);
        
    % Project the new sample
    test = bsxfun(@minus,sp(idxTest,:),meanTrain);
    allS(idxTest,:) = test * ll;

    % What about a knn predictive model?
    mdl = fitcknn(ss,grpIdx,'NumNeighbors',3,'Standardize',0);
    pre = predict(mdl,allS(idxTest,:));
    allP(idxTest,:) = [allIdx(idxTest),pre];

    % Add to the confusion matrix...
    fx = find(idxTest);
    for r = 1:numel(pre)        
        cm(allIdx(fx(r)),pre(r)) = cm(allIdx(fx(r)),pre(r)) + 1;        
    end
    
    % Draw...
    if verb
        plotScores(f1,ss,grpIdx,allS,numG,idxTest);
    end

    waitbar(n/k,wb);
end

% Rescale the loadings by taking the mean
allL = allL / numS;

delete(wb);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotScores(f1,ss,grpIdx,allS,numG,idxTest)

% Get and draw the confidence ellipses...
%clear ci
for r = 1:numG
    [ci(r).ell,ci(r).xy,ci(r).ab] = error_ellipse(ss(grpIdx == r,1:2),95);
end

figure(f1);
cla reset;

hold on;
scatter(ss(:,1),ss(:,2),80,grpIdx,'o','filled',...
    'MarkerEdgeColor',[0.8 0.8 0.8]);

for r = 1:numG
    plot(ci(r).ell(:,1),ci(r).ell(:,2),'Color','k','LineWidth',1);
    scatter(ci(r).xy(1),ci(r).xy(2),50,'rx');
end

% Formattings
box on;
set(gca,'FontSize',14);
xlabel('LV1','FontSize',16,'FontWeight','bold');
ylabel('LV2','FontSize',16,'FontWeight','bold');

scatter(allS(idxTest,1),allS(idxTest,2),100,'r','o','filled',...
    'MarkerEdgeColor','k');        

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%