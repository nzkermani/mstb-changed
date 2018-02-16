function falsePositives
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

numTests = 100;
numPoints = 1000;
sampleSize = 0.1;

% Generate data
rr = randn(numPoints,1);

% Store pvalues
pp = NaN(numTests,1);

% Store significant results
signif = zeros(numPoints * sampleSize,numTests * 0.1);
i = 0;

for n = 1:numTests
    
    idx = randperm(1000,numPoints*sampleSize);
    
    [~,pp(n,1)] = ttest(rr(idx));
    
    if pp(n,1) < 0.05
        i = i + 1;
        signif(:,i) = rr(idx);
    end
    
end

% Trim out the crap
signif = signif(:,1:i);

figure;
subplot(1,2,1); 
hist(rr,100);

subplot(1,2,2);
boxplot(signif);


end

