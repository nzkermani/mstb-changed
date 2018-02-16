% Make some plots for the paper

% Align the two images
opts.ppmTol = 5;
opts.handBand = eval(['@(z) (' num2str(opts.ppmTol) '* z / 1e6)']);            
opts.handGap  = eval(['@(z) (' num2str(opts.ppmTol) '* z(:,1) / 1e6)']);
opts.handDist = @(R,S) abs(sum((R-S),2));

sp1 = [dpn.d1.mz' ones(numel(dpn.d1.mz),1)];
sp2 = [dpn.d2.mz' ones(numel(dpn.d2.mz),1)];

[j,k] = samplealign2(...
    sp1,...
    sp2,...
    'Band',opts.handBand,...
    'Gap',opts.handGap,...
    'Distance',opts.handDist,...
    'Quantile',[],...
    'SHOWCONSTRAINTS',false,...
    'SHOWNETWORK',false,...
    'SHOWALIGNMENT',false);

% Make new
mz1 = dpn.d1.mz(j);
mz2 = dpn.d2.mz(k);
sp1 = dpn.d1.sp(:,:,j);
sp2 = dpn.d2.sp(:,:,k);

sp1 = bsxfun(@rdivide,sp1,nansum(sp1,3));
sp2 = bsxfun(@rdivide,sp2,nansum(sp2,3));

% Now determine correlation / deviation of each image
crr = zeros(numel(mz1),1);
dff = zeros(numel(mz1),2);
for n = 1:numel(mz1)
    
    t1 = sp1(:,:,n);
    t2 = sp2(:,:,n);
    t1 = t1(:);
    t2 = t2(:);
    
    fx = t1 > 0 & t2 > 0;
    
    crr(n,1) = corr(t1(fx),t2(fx));
    
    tmp = t1 - t2;
    dff(n,:) = [mean(tmp(fx)) std(tmp(fx))];

end


ld1 = dpn.fuse.mva.PCA.ld.d1(j,1);
ld2 = dpn.fuse.mva.PCA.ld.d2(k,1);

figure; hold on;
stem(mz1,ld1,'MarkerSize',0.01,'LineWidth',2);
stem(mz1,ld2,'MarkerSize',0.01,'LineWidth',2);

figure;
fg = 100*(ld1-ld2)./ld1;
mask = fg < 1e4;
stem(mz1(mask),fg(mask));
