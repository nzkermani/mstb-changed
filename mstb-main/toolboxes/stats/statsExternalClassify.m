function statsExternalClassify(~,~,fig,window)
% statsExternalClassify - load the samples and classify against the model

% Get the name of the test set file
test = get(window.load,'UserData');
if isempty(test)
    disp('No test set data found');
    return
end

% Get the guidata
sts = guidata(fig.fig);

% Determine that a model has been built...
model = sts.res.mmc;

% Open the test set file and determine the number of files
tmp = open([test.path test.file]);
numF = size(tmp.res.data,2);

% Which labels are the 'true' values? These were chosen in classCompare
val = get(window.classCompare,'Value');
str = get(window.classCompare,'String');
compMeta = tmp.res.meta.(str{val});

% Correct 'B1 ' and 'B1' being different
fx = strcmp(compMeta,'B1 ');
compMeta(fx) = {'B1'};

% Somewhere to save the projected scores, that way we can predict all
% together once the scores have been determined
predScores = zeros(numF,size(model.ss,2));

% Calculate the projected scores
for n = 1:numF
    predScores(n,:) = externalProject(sts.proc,sts.res.mmc,tmp.res.data(n));
end

% Now here we need to classify the samples using one of the available
% methods. This is easily done for centroid mode and knn. Note that we want
% to be able to classify using the appropriate number of dimensions
%[pre] = externalClassify(sts.res.mmc,predScores,'knn',3);
[pre] = externalClassify(sts.res.mmc,predScores,'centroid',3);

% Generate a confusion matrix accordingly - note that the labels are not
% designed to match perfectly for actual and predicted class...
[unqA,~,indA] = unique(compMeta);   % actual class of samples
[unqP,~,indP] = unique(pre);        % predicted class
cm = zeros(numel(unqA),numel(unqP));
for n = 1:size(pre,1)
    
    cm(indA(n),indP(n)) = cm(indA(n),indP(n)) + 1;
    
end

makeConfMatPredict(cm,unqA,unqP);


figure; hold on;

% MODEL
[unqM,~,indM] = unique(sts.res.mmc.grp);
lgX = zeros(numel(unqM)+numel(unqA),1);
cols = parula(numel(unqM));
for n = 1:numel(unqM)
    fx = indM == n;
    lgX(n) = scatter(sts.res.mmc.ss(fx,1),sts.res.mmc.ss(fx,2),...
        120,cols(n,:),'o','filled',...
        'MarkerEdgeColor','k');
end

% Misclassified ones
fx = indA ~= indP;
scatter(predScores(fx,1),predScores(fx,2),100,'o','red')

lg2 = zeros(numel(unqA),1);
cols = parula(numel(unqA));
for n = 1:numel(unqA)
    fx = indA == n;
    lgX(n+numel(unqM)) = scatter(predScores(fx,1),predScores(fx,2),...
        50,cols(n,:),'d','filled',...
        'MarkerEdgeColor','k');
end



% Place legend
legend(lgX,[unqM; unqA])

box on;

set(gca,'FontSize',14);
xlabel('LV1','FontSize',16,'FontWeight','bold');
ylabel('LV2','FontSize',16,'FontWeight','bold');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [testScores] = externalProject(train,model,test)
% This function needs to perform alignment of the spectra to match the
% training set, and then project onto the model accordingly...

% Align test file to mz/sp from training set
res = 0.1;
[~,data] = dbBinningFixed(test,[min(train.var.mz) max(train.var.mz)],res);
testSp = data.al;

% Now we need to process testSp as was done for the model

% Normalise (we'll cheat with TIC)
scFac = mean(nansum(train.sp,2));
testSp = scFac * testSp / sum(scFac);

% Transform
if train.opts.doLog
    testSp = log(testSp + train.opts.logOS);
end

% Scaling not currently supported!
if strcmp(train.opts.scale,'UV')
    testSp = zscore(testSp);
end

% Now we need to project onto the model I guess?
meanX = nanmean(train.sp,1);
testSp = bsxfun(@minus,testSp,meanX);
testScores = testSp * model.ll;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cl] = externalClassify(model,predSc,method,k)
% Classify the samples!
        

% Determine unique classes and indices
[unq,~,ind] = unique(model.grp);
numG = numel(unq);

% How many components to use? This is typically fixed at numG-1
numDV = numG - 1;

switch method
    
    case 'knn'

        % Fit kNN model
        mdl = fitcknn(model.ss(:,1:numDV),ind,'NumNeighbors',k,'Standardize',0);
        
        % Now predict!
        pre = predict(mdl,predSc(:,1:numDV));
        
        % Put into a cell array of classes

    case 'centroid'
        
        
        % Determine centroids of each class in model
        cntr = zeros(numG,numDV);
        for n = 1:numG
            fx = ind == n;
            cntr(n,:) = nanmean(model.ss(fx,1:numDV),1);            
        end
        
        % Distance to centroids
        tmpDist = pdist2(predSc(:,1:numDV),cntr);
        
        % Closest one wins
        [~,pre] = min(tmpDist,[],2);
            
        
    otherwise
        
end

% Convert the numeric entries from 'pre' to text ones corresponding to
% those in model.grp
cl = repmat({'NA'},[numel(pre) 1]);
for n = 1:numG
    
    fx = pre == n;
    cl(fx) = unq(n);
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%