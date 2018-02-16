function [ lc ] = lpolcMMC(sp,grp,cvp)
%lpolcMMC - leave group (patient) out MMC with learning curve function. 
%Modified version to increase training size. This is a supervised analysis
%method that runs MMC with each group specified in grp being projected
%entirely at the same time. This reduces bias and is more analogous to what
%would happen under in vivo conditions.

% Define parameters for the learning curve
minPC = 2; % start point, spectra per class
stepC = 1; % step size, i.e. increase minPC by this each iteration
numRepl = 30; % repeat

% Create a waitbar
wb = waitbar(0,'kFold CV');

% Create the indices for which groups to be omitted in turn
[cvUnq,~,cvIdx] = unique(cvp);
numCV = numel(cvUnq);

% How many DVs to calculate?
[unq,~,allIdx] = unique(grp);
numG = numel(unq);
numDV = max([2 numG-1]);

% Need to define the rough step sizes based on the number of spectra...
clSizes = hist(allIdx,1:numG);
lcSize = minPC:stepC:min(clSizes);
numStep = numel(lcSize);
    
% Store the classification results
clCorr = zeros(numStep,numRepl,numCV);
clInco = zeros(numStep,numRepl,numCV);
clNorm = zeros(numCV,2);

% Loop through the folds
for n = 1:numCV
    
    % Define the training/testing set indices
    idxTrain = cvIdx ~= n;
    idxTest  = cvIdx == n;
    
    % Need to separate this for subsequent indexing
    spTrain = sp(idxTrain,:);
    
    % Define the groupings within the training set
    [~,~,grpIdx] = unique(grp(idxTrain));
    
    % Here I want to perform a full-sized classification
    corr = runMMC(spTrain,grpIdx,numDV,sp(idxTest,:),allIdx(idxTest,:));
    clNorm(n,:) = [sum(corr == 1) sum(corr == 0)];
        
    % We need to split up the training set into further chunks
    for r = 1:numStep
       
        % Need to perform for a certain number of replicates to ensure that
        % that the training sets are as random as possible
        for s = 1:numRepl
        
            % Now here we can prepare the subset training set
            [lcIdx] = generateSubset(grpIdx,lcSize(r),numG);
        
            % Here we run the function...
            [corr] = runMMC(spTrain(lcIdx,:),...% training data
                grpIdx(lcIdx,:),...% training class
                numDV,...
                sp(idxTest,:),...% test data
                allIdx(idxTest,:));% test class
            
            % Save the results...
            clCorr(r,s,n) = sum(corr == 1);
            clInco(r,s,n) = sum(corr == 0);
                
        end
        
    end
    
    waitbar(n/numCV,wb,cvUnq{n});
end

delete(wb);

% Calculate the important values and such like from this... We take the
% mean value over each iteration and then the mean over all omitted
% groups...
coCl = sum(clCorr,3);
inCl = sum(clInco,3);

vals = 100 * coCl ./ (coCl + inCl);
per = mean(vals,2);
dev = std(vals,[],2);

best = sum(clNorm(:,1)) / sum(clNorm(:)) * 100;

% Save the learning curve outputs...
lc.step = lcSize;
lc.class = per;
lc.dev = dev;
lc.best = best;

figure;
errorbar(lcSize,per,dev,...
    '-s','MarkerSize',10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor','red',...
    'LineWidth',2);

% hold on;
% line([0 lcSize(end)+stepC],[best best],...
%     'LineWidth',2,...
%     'Color','k',...
%     'LineStyle','--');

set(gca,'FontSize',14);
xlabel('Training Set Size (per class)','FontSize',18,'FontWeight','bold');
ylabel('Classification performance / % ±1\sigma','FontSize',18,'FontWeight','bold');

xlim([0 lcSize(end)+stepC]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lcIdx] = generateSubset(grpIdx,sampleSize,numG)
% Randomly partition the training set...

% Store indices in here...
lcIdx = zeros(sampleSize,numG);

% For each group
for n = 1:numG
    
    % Find indices of group n
    fx = find(grpIdx == n);
    
    % Get the random indices
    fy = randperm(numel(fx),min([sampleSize numel(fx)]));
    
    % Save to the thing
    lcIdx(1:numel(fy),n) = fx(fy);
end

% Make a single column
lcIdx = lcIdx(:);

% Remove any zero values from here
mask = lcIdx ~= 0;
lcIdx = lcIdx(mask);

end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [corr] = runMMC(sp,grpIdx,numDV,testSp,testClass)

% Run the MMC on the training set
[~,ss,~,ll] = recursiveMmcLda(sp,grpIdx,numDV);

% Determine the mean spectrum of the training set
meanTrain = nanmean(sp,1);

% Project the new samples to get scores
testSp = bsxfun(@minus,testSp,meanTrain);
scTest = testSp * ll;

% What about a knn predictive model?
mdl = fitcknn(ss,grpIdx,'NumNeighbors',3,'Standardize',0);
pre = predict(mdl,scTest);

% Vector of correctness
corr = testClass == pre;

% figure; hold on;
% scatter(ss(:,1),ss(:,2),80,grpIdx,'o','filled');
% scatter(scTest(:,1),scTest(:,2),80,'k','d','filled');


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%