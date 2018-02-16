function [mean1,mean2] = dpnPaperLoadings(dpn)
%dpnPaperLoadings - display big loadings for table / publication

% These are the loadings...
ld = dpn.fuse.mva.MMC.ld;

figure; hold on;
stem(ld.mz1,ld.d1(:,1),'MarkerSize',0.01,'LineWidth',2,'Color','b');
stem(ld.mz2,ld.d2(:,1),'MarkerSize',0.01,'LineWidth',2,'Color','r');

% Find the biggest ones from each...
[~,ix1] = sort(abs(ld.d1(:,1)),'descend');
[~,ix2] = sort(abs(ld.d2(:,1)),'descend');

% How many do we want?
maxL = 10;

% Extract mz values
[sort(ld.mz1(ix1(1:maxL)))' sort(ld.mz2(ix2(1:maxL)))']
[(ld.mz1(ix1(1:maxL)))' (ld.mz2(ix2(1:maxL)))']

% Can we determine the fold change?
[a,b,~,~] = dpnAnnotationExtract(dpn);
a(a == 15) = 0;
mask = a ~= 0;
histID = b(mask);
[a,~,c] = unique(histID);
a

% Reshape original data...
sz1 = size(dpn.d1.sp);
sp1 = reshape(dpn.d1.sp,[sz1(1)*sz1(2) sz1(3)]);
sz2 = size(dpn.d2.sp);
sp2 = reshape(dpn.d2.sp,[sz2(1)*sz2(2) sz2(3)]);

% Extract useful pixels...
sp1 = sp1(mask,:);
sp2 = sp2(mask,:);

% Now we TIC normalise only (no log transform)
sp1 = bsxfun(@rdivide,sp1,nansum(sp1,2)) * 1000;
sp2 = bsxfun(@rdivide,sp2,nansum(sp2,2)) * 1000;

mean1 = zeros(numel(a),size(sp1,2));
mean2 = zeros(numel(a),size(sp2,2));
for n = 1:numel(a)
    mean1(n,:) = nanmean(sp1(c == n,:),1);
    mean2(n,:) = nanmean(sp2(c == n,:),1);
end

%mean1 = [nanmean(sp1(c == 1,:),1); nanmean(sp1(c == 2,:),1)];
%mean2 = [nanmean(sp2(c == 1,:),1); nanmean(sp2(c == 2,:),1)];

return

% Annotation
[lm,ass] = annotateMZ(dpn.d1.mz,...
    'Polarity','pos',...
    'Adduct',{'M+H','M+Na','M+K','M+NH4'},...
    'Tolerance',10,...
    'Database','dipa');
annotateOP2(dpn.d1.mz,sp1,histID,ass,lm,...
    '/Users/jmckenzi/Desktop/Anno-Positive.txt',0.05);

[lm,ass] = annotateMZ(dpn.d2.mz,...
    'Polarity','neg',...
    'Adduct',{'M-H'},...
    'Tolerance',10,...
    'Database','dipa');
annotateOP2(dpn.d2.mz,sp2,histID,ass,lm,...
    '/Users/jmckenzi/Desktop/Anno-Negative.txt',0.05);

end

