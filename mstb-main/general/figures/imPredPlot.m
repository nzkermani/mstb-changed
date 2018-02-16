function imPredPlot( train,test )
%imPredPlot - a follow on from imagePrediction.m
%
% Various plots can be generated from this, possibly



% Plot the MMC scores of the annotated pixels
plotAnno(train)




end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotAnno(train)
% A function to plot the annotations of the image, and the scores and the
% images for classification

% Reshape the scores
sc = reshape(train.score,[train.sz(1)*train.sz(2) size(train.score,3)]);

% Determine non-zero annotated indices
anno = bsxfun(@times,train.anno,1:size(train.anno,2));
anno = nansum(anno,2);

% Which indices were annotated?
fx = anno > 0 & anno ~= 4;

figure; hold on;
scatter(sc(fx,1),sc(fx,2),80,anno(fx),'o','filled');

figure; hold on;
scatter(sc(fx,1),sc(fx,3),80,anno(fx),'o','fiiled');


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%