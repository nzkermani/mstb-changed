function [ allS,cm,allP,unq,allL ] = looMMC(sp,grp)
%looMMC - leave one pixel out cross validated MMC
%
% For the spectral data and its groupings, we need to iteratively perform
% leave one pixel out MMC and then predict the left out pixel. Perhaps by
% default we will use 3-nearest neighbours. Also try to present the
% confusion matrix

verb = false;
if verb
    f1 = figure;
    ax = axes('Parent',f1);
end

% Create a waitbar
wb = waitbar(0,'LOO CV');

% Number of spectra
[numS,numV] = size(sp);
trainVec = (1:numS)';

% How many DVs to calculate?
[unq,~,allIdx] = unique(grp);
numG = numel(unq);
numDV = max([2 numG-1]);

% Somewhere to store the predicted scores...
allS = zeros(numS,numDV);
allP = zeros(numS,2);
allL = zeros(numV,numDV);

% Confusion matrix...
cm = zeros(numG,numG);

% Run a simple, single analysis using ALL available observations
[~,~,~,totLL] = recursiveMmcLda(sp,allIdx,numDV);

% Loop
for n = 1:numS
    
    % Define the training set
    idx = trainVec ~= n;
    
    % Define the groupings
    [~,~,grpIdx] = unique(grp(idx));
    
    % Run the MMC on the training set
    [~,ss,~,ll] = recursiveMmcLda(sp(idx,:),grpIdx,numDV);
    
    % We should mirror the scores in the case of 2 groups in order that
    % they are compatible.  This prevents the scores from looking stupid if
    % the first axes gets swapped around. I'm not sure how this will be
    % shown if there are more than 2 groups in the analysis
    [~,vals] = alignCVdirsDR(totLL,ll);
    ll = bsxfun(@times,ll,vals);
    ss = bsxfun(@times,ss,vals);
            
    
    % Save the loadings
    allL = allL + ll;
    
    % Determine the mean spectrum of the training set
    meanTrain = nanmean(sp(idx,:),1);
    
    % Project the new sample
    test = bsxfun(@minus,sp(~idx,:),meanTrain);
    allS(n,:) = test * ll;
    
    % Classify the sample...
    
    % What about a knn predictive model?
    mdl = fitcknn(ss(:,1:numG-1),grpIdx,'NumNeighbors',3,'Standardize',0);
    pre = predict(mdl,allS(n,1:numG-1));
    allP(n,:) = [allIdx(n),pre];
    
    % Add to the confusion matrix...
    cm(allIdx(n),pre) = cm(allIdx(n),pre) + 1;
    
    if verb
        
        % Get and draw the confidence ellipses...
        clear ci
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
        
        scatter(allS(n,1),allS(n,2),200,[1 1 1],'d','filled',...
            'MarkerEdgeColor','k');        
        scatter(allS(n,1),allS(n,2),80,allIdx(n),'d','filled',...
            'MarkerEdgeColor','k');
        
        % Line plot
        for r = 1:numG            
            plot([allS(n,1) ci(r).xy(1)],[allS(n,2) ci(r).xy(2)],'--',...
                'Color','k');
        end
        
        

        %pause(2);
    end
    
    waitbar(n/numS,wb);
end

% Rescale the loadings by taking the mean
allL = allL / numS;

delete(wb);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%