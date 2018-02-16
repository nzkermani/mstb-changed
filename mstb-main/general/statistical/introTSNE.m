function [ mapTrain ] = introTSNE(xy,idx)
%introTSNE - starting point for tSNE analysis

% Determine the input arguments
[xy,idx,train] = generateData;
% xy = rand(2000,10);

% Set some tSNE parameters
numDims = 2;
initDims = 2;
perplex = 100;
%[xy,idx,train] = generateData;
%xy = rand(2000,10);

% Set some tSNE parameters
numDims = 4;
initDims = 20;
perplex = 30;
theta = 0.5;

% Train the tSNE network
tic
mapTrain = tsne(xy,[],numDims,initDims,perplex);
toc

if nargin == 1
    return
end    

% Plot test embedding
figure('Units','normalized','Position',[0.1 0.2 0.8 0.3]);
% ax(1) = subplot(1,2,1);
% scatter3(xy(:,1),xy(:,2),xy(:,3),30,idx,'o','filled');
% 
% ax(2) = subplot(1,2,2);
scatter(mapTrain(:,1),mapTrain(:,2),80,idx,'o','filled',...
    'MarkerEdgeColor','k');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xy,idx,tr] = generateData
% Generate random data for testing...

clustCent = [-10 -10 10 2; -5 -10 -10 0.5; 5 10 0 1; 10 5 5 1; 0 0 0 1];
numClust = size(clustCent,1);
sizeClust = 50;
numDim = size(clustCent,2) - 1;
xy = cell(numClust,3);
for n = 1:numClust
    
    tmp = randn(sizeClust,numDim) * clustCent(n,end);
    xy{n,1} = bsxfun(@plus,tmp,clustCent(n,1:end-1));
    xy{n,2} = ones(sizeClust,1) * n;
    
    % Determine test/train from each cluster...
    xy{n,3} = rand(sizeClust,1) > 0.8;
    
    

end

idx = vertcat(xy{:,2});
tr = vertcat(xy{:,3}) == 0;
xy = vertcat(xy{:,1});

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%