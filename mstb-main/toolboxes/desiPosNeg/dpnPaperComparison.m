function [ op ] = dpnPaperComparison(dpn,verb)
%dpnPaperComparison - analyses of A30 sample

% Draw figures unless explicitly not required
if nargin == 1
    verb = true;
end
if ~islogical(verb)
    verb = true;
end

% Let's just plot the variables first...
s1 = size(dpn.d1.sp);
s2 = size(dpn.d2.sp);

sp1 = reshape(dpn.d1.sp,[s1(1)*s1(2) s1(3)]);
sp2 = reshape(dpn.d2.sp,[s2(1)*s2(2) s2(3)]);

m1 = nanmean(sp1,1);
m2 = nanmean(sp2,1);

% Now let's try to match the smallest to the first...
pks = NaN(s2(3),5);
[~,srtIdx] = sort(m1,'descend');
mzVec = dpn.d2.mz;
for i = 1:s1(3)
    
    n = srtIdx(i);
    
    % Find the nearest available peak...
    [ds,fx] = min(abs(mzVec - dpn.d1.mz(n)));
    
    if ds < 0.05
        pks(fx,:) = [fx n mzVec(fx) dpn.d1.mz(n) mzVec(fx)-dpn.d1.mz(n)];
        mzVec(fx) = NaN;
    else
        % Do nothing
    end
    
    
end

fx = isnan(pks(:,1));
fx = find(fx);

fy = isnan(pks(:,2));
%fy = find(fy);

% Draw a figure to show peak matching
if verb
    figure; hold on;
    stem(dpn.d1.mz,m1,'b');
    stem(dpn.d2.mz,-m2,'r');
    scatter(dpn.d2.mz(fx),m2(fx),80,'k','o','filled')
end

% New matrices corresponding to matched variables
n1 = sp1(:,pks(~fy,2));
n2 = sp2(:,pks(~fy,1));
crr = NaN(size(n1,1),1);
for n = 1:size(n1,1)
    crr(n,1) = corr(n1(n,:)',n2(n,:)');
end
crr = reshape(crr,[s1(1) s1(2)]);

% Correlation plot figure
if verb
    figure; 
    ax = imagesc(crr);
    axis square
    axis off
    cmap = hot(100);
    colormap(cmap);

    cb = colorbar;
    set(gca,'FontSize',16);
    ylabel(cb,'Correlation Coefficient','FontSize',16,'FontWeight','bold');
end


% Want an output...
op.crr = crr;
op.n1 = reshape(n1,[s1(1) s1(2) size(n1,2)]);
op.n2 = reshape(n2,[s2(1) s2(2) size(n2,2)]);
op.mz1 = dpn.d1.mz(pks(~fy,2));
op.mz2 = dpn.d2.mz(pks(~fy,1));
op.interp1 = dpn.d1.isInterp;
op.interp2 = dpn.d2.isInterp;

end

