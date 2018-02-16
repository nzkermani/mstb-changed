function interpCompare(orbData,mz1,op1,mz2,op2)
% Need to compare a single spectrum from the orbitrap data, i.e. raw data
% to the interpolated forms (both Da and ppm style)...

mzRange = [885.4 885.6];



raw = orbData{1};
numS = size(raw,2);

% Spectral masks
mask1 = mz1 > mzRange(1) & mz1 < mzRange(2);
mask2 = mz2 > mzRange(1) & mz2 < mzRange(2);

normFac = [max(op1(1,mask1)) max(op2(1,mask2))]



figure; hold on;



% Have to loop through and plot them all
for n = 1:numS
    
    mask = raw{n}(:,1) > min(mzRange) & raw{n}(:,1) < max(mzRange);
    
    x = raw{n}(mask,1);
    y = mean(normFac) * raw{n}(mask,2) / max(raw{n}(mask,2));
    
    
    
    plot(x,y,'b');
    
end

% Add in the interpolated forms
stem(mz1(mask1),op1(1,mask1),'r',...
    'MarkerFaceColor','red',...
    'MarkerEdgeColor','black',...
    'MarkerSize',10);

stem(mz2(mask2),op2(1,mask2),'g',...
    'MarkerFaceColor','green',...
    'MarkerEdgeColor','black',...
    'MarkerSize',10);


end