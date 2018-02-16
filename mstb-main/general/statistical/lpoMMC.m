function [ allS,cm,allP,unq,allL,roc ] = lpoMMC(sp,grp,cvp)
%lpoMMC - leave group (patient) out MMC. THis is a supervised analysis
%method that runs MMC with each group specified in grp being projected
%entirely at the same time. This reduces bias and is more analogous to what
%would happen under in vivo conditions.

% Define the classification method
classMethod = 'centroid';

verb = false;
if verb
    f1 = figure;
    ax = axes('Parent',f1);
end

% Create a waitbar
wb = waitbar(0,'kFold CV');

% Create the indices for which groups to be omitted in turn
try
    [cvUnq,~,cvIdx] = unique(cvp);
catch
    cvp = cell2mat(cvp);
    [cvUnq,~,cvIdx] = unique(cvp);
end
numCV = numel(cvUnq);

% Number of spectra
[numS,numV] = size(sp);

% How many DVs to calculate?
[unq,~,allIdx] = unique(grp);
numG = numel(unq);
numDV = max([2 numG-1]);
calcDV = numG - 1;

% Somewhere to store the predicted scores...
allS = zeros(numS,numDV);
allP = zeros(numS,2);
allL = zeros(numV,numDV);

% Run a simple, single analysis using ALL available observations
[~,totSS,~,totLL] = recursiveMmcLda(sp,allIdx,numDV);

% Confusion matrix...
cm = zeros(numG,numG);

% Loop through the folds
for n = 1:numCV
    
    % Define the training/testing set indices
    idxTrain = cvIdx ~= n;
    idxTest  = cvIdx == n;
    
    % Define the groupings
    [grpUnq,~,grpIdx] = unique(grp(idxTrain));
    numUnq = numel(grpUnq);
    
    % Run the MMC on the training set
    [~,ss,~,ll] = recursiveMmcLda(sp(idxTrain,:),grpIdx,numDV);
    
    % This function (replacing the previous one) needs to map the
    % calculated scores on to the overall scores to ensure consistency.
    % This is because axes can be flipped which will ruin the average
    % scores to be presented at the end of the analysis...
    [~,vals] = alignCVdirsDR(totLL,ll);
    ll = bsxfun(@times,ll,vals);
    ss = bsxfun(@times,ss,vals);
    
    % Save the loadings
    allL = allL + ll;
    
    % Determine the mean spectrum of the training set
    meanTrain = nanmean(sp(idxTrain,:),1);
        
    % Project the new samples
    test = bsxfun(@minus,sp(idxTest,:),meanTrain);
    allS(idxTest,:) = test * ll;
    

    % How do we wish to classify the samples?
    switch classMethod
        
        case 'centroid'

            % Simple centroid distances
            mhld = zeros(numUnq,sum(idxTest));
            for r = 1:numUnq
                mx = grpIdx == r;

                % Determine centroid
                tmpCent = nanmean(ss(mx,1:calcDV),1);

                % Determine distance between centroid and predicted observations
                tmpDist = pdist2(allS(idxTest,1:calcDV),tmpCent);

                % Save to mhld
                mhld(r,:) = tmpDist';

            end

            % Find minimum distance
            [~,pre] = min(mhld,[],1);
            allP(idxTest,:) = [allIdx(idxTest),pre'];

        case 'mahal'
            
            % Use the Matlab provided mahalanobis distance
            pre = classify(allS(idxTest,:),ss,grpIdx,'mahalanobis');
            allP(idxTest,:) = [allIdx(idxTest),pre];
    
        case 'knn'

            % What about a knn predictive model?
            mdl = fitcknn(ss,grpIdx,'NumNeighbors',3,'Standardize',0);
            pre = predict(mdl,allS(idxTest,:));
            allP(idxTest,:) = [allIdx(idxTest),pre];

    end

    % Add to the confusion matrix...
    fx = find(idxTest);
    for r = 1:numel(pre)        
        cm(allIdx(fx(r)),pre(r)) = cm(allIdx(fx(r)),pre(r)) + 1;        
    end
    
    % Draw...
    if verb
        plotScores(f1,ss,grpIdx,allS,numG,idxTest);
        title(cvUnq{n});
    end

    if isnumeric(cvUnq)
        waitbar(n/numCV,wb,cvUnq(n));
    else
        waitbar(n/numCV,wb,cvUnq{n});
    end
end

% If there are two curves, then we do the ROC curve for the non-CV and CV
% scores
if numG == 2
    
    roc.x = [0:0.01:1]';
    [roc.xNoCV,roc.yNoCV,~,roc.aucNoCV] = perfcurve(allIdx,totSS(:,1),2,'XVals',roc.x);
    [roc.xCV,roc.yCV,~,roc.aucCV] = perfcurve(allIdx,allS(:,1), 2,'XVals',roc.x);
   
else
    roc = [];
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