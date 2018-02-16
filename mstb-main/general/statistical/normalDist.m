function normalDist
%normalDist - play functions for plots and stuff...

numD = 5;

sigma = 1:1:numD;
mu = linspace(1,100,numD);

% Determine histograms...
cols = parula(numD);
figure; hold on;

for n = 1:numD    
    
    rn = normrnd(mu(n),sigma(n),500,1);
    
    [a,b] = hist(rn,-10:(max(mu)*1.2));
    
    hh = bar(b,a,'FaceColor',cols(n,:));
    %set(hh.FaceColor,cols(n,:));
    
end


end

