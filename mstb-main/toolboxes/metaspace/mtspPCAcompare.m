function [ output_args ] = mtspPCAcompare(res,dpn)
%mtspPCAcompare - two plots for PCA...


dpn

% Extract annotated regions for this file
[grp,histID,~] = desiAnnotationExtract(dpn,false);
mask = grp > 0;

% Normalise the data somehow
doLog = true;
[sp1,flag1] = prepare([],dpn.mtspOLD,doLog);
[sp2,flag2] = prepare([],dpn.mtsp,doLog);    

% Run PCA
[ll1,ss1,~] = pca(sp1,'Economy',true);
[ll2,ss2,~] = pca(sp2,'Economy',true);

% Consider reshaping if in image form
if flag1
    ss1 = reshape(ss1(:,1:10),[size(dpn.mtsp.sp,1) size(dpn.mtsp.sp,2) 10]);
    ss2 = reshape(ss2(:,1:10),[size(dpn.mtsp.sp,1) size(dpn.mtsp.sp,2) 10]);
end

% Determine the two categories of variables:
% 1 - missing from new file
% 2 - new in the new file
[tab] = mtspCSVcompare(res,dpn);
fdrs = cell2mat(tab(:,7:8));
fdrs(isnan(fdrs)) = Inf;

% Old
fx = fdrs(:,1) <= 0.1;
fy = fdrs(:,2) > 0.1;
fOO = fx & fy;

% New
fx = fdrs(:,2) <= 0.1;
fy = fdrs(:,1) > 0.1;
fNN = fx & fy;

[mOO] = match2mass(tab(fOO,:),dpn.mtspOLD);
[mNN] = match2mass(tab(fNN,:),dpn.mtsp);

% Figure - two axes to show the scores adjacent to each other
figure;

if flag1
    subplot(2,2,1); imagesc(imScale(ss1(:,:,1:3)));
    axis off;
    subplot(2,2,2); imagesc(imScale(ss2(:,:,1:3)));
    axis off;
    
    mzo = mtspForm2MZ(dpn.mtspOLD.annos);
    subplot(2,2,3); hold on;
    stem(mzo,ll1(:,1),'k','MarkerSize',0.01);
    stem(mzo(mOO),ll1(mOO,1),'r','MarkerFaceColor','r');
    box on;
    set(gca,'FontSize',14);
    xlabel('m/z','FontSize',16);
    ylabel('PC1 Loading','FontSize',16);
    
    mzo = mtspForm2MZ(dpn.mtsp.annos);
    subplot(2,2,4); hold on;
    stem(mzo,ll2(:,1),'k','MarkerSize',0.01);
    stem(mzo(mNN),ll2(mNN,1),'r','MarkerFaceColor','r');
    box on;
    set(gca,'FontSize',14);
    xlabel('m/z','FontSize',16);
    ylabel('PC1 Loading','FontSize',16);
    
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [match] = match2mass(newtab,data)
% Match molecular formulae in one list to those in the dpn structure

match = false(numel(data.mz),1);

% Can we match the data...
for n = 1:size(newtab,1)
    
    thisForm = newtab{n,1};
    thisAdct = newtab{n,2};
    
    % Match to existing formulae
    fx = strcmp(data.annos(:,1),thisForm);
    fy = strcmp(data.annos(:,2),thisAdct);
    fz = fx & fy;
    if sum(fz) == 1
        i = find(fz);
        match(i,1) = true;
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sp,flag] = prepare(mask,data,doLog)
% Prepare the data...

% Reshape
sp = data.sp;
sz = size(sp);
sp = reshape(sp,[sz(1)*sz(2) sz(3)]);

if ~isempty(mask)
    sp = sp(mask,:);
    flag = false;
else
    flag = true;
end

% Normalise
[sp] = jsmNormalise(sp,'pqn-median',0,0,[]);

% Transform
if doLog
    logOS = nanmedian(sp(sp > 0));
    minOS = min(sp(:));
    if minOS + logOS < 1
        logOS = 1;
    end
    sp = log(sp + logOS);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

