function [ output_args ] = twc(aa,bb,cc)
%twc - three way comparison of data that has been differently interpolated
%
% Do PCA and then find a good way to present the results

grp = [1 1 1 1 1 ...
    2 2 2 2 2 ...
    3 3 3 3 3 ...
    4 4 4 4 4 ...
    5 5 5 5 5 ...
    6 6 6 6 6 ...
    7 7 7 7 7 ...
    8 8 8 8 8 ...
    9 9 9 9 9 ...
    10 10 10 10 10 ...
    11 11 11 11 11 ...
    12 12 12 12 12];

% NB sample 54 is poorly aligned. Just ditch it
numO = size(aa.op,1);
idx = true(numO,1);
idx(54) = false;

% We perform TIC norm beforehand on all spectr
aa.op = bsxfun(@rdivide,aa.op,nansum(aa.op,2));
bb.op = bsxfun(@rdivide,bb.op,nansum(bb.op,2));
cc.op = bsxfun(@rdivide,cc.op,nansum(cc.op,2));

% Now do the PCA
[laa,saa,~] = princomp(aa.op(idx,:),'econ');
[lbb,sbb,~] = princomp(bb.op(idx,:),'econ');
[lcc,scc,~] = princomp(cc.op(idx,:),'econ');

% Do procrustes analysis, using scc as a fixed point
[~,paa] = procrustes(scc,saa);
[~,pbb] = procrustes(scc,sbb);


% Loadings plot for scumsum
cu1 = cumsum(abs(laa(:,1)));
cu2 = cumsum(abs(lbb(:,1)));
cu3 = cumsum(abs(lcc(:,1)));

figure; hold on;
plot(aa.mz,cu1/max(cu1),'-ro');
plot(bb.mz,cu2/max(cu2),'-go');
plot(cc.mz,cu3/max(cu3),'-bo');




cols = eye(3);
markers = {'o','s','d'};

figure; hold on;

scatter(paa(:,1),paa(:,2),60,grp(idx),markers{1},'filled');
scatter(pbb(:,1),pbb(:,2),60,grp(idx),markers{2},'filled');
scatter(scc(:,1),scc(:,2),60,grp(idx),markers{3},'filled');







end

