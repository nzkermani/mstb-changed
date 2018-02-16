function [ output_args ] = kendrickClass(mz,cl)
%kendrickClass - visualisation plot

% Calculate the KMD
[nom,kmd] = kendrick(mz,[],'ch2');

% Now scatter through each type...
[unq,~,ind] = unique(cl);
numC = numel(unq);

cols = parula(numC);
ax = zeros(numC,1);

figure; hold on;

for n = 1:numC
    
    
    fx = ind == n;
    
    ax(n) = scatter(nom(fx),kmd(fx),90,cols(n,:),'o','filled',...
        'MarkerEdgeColor',[0.8 0.8 0.8]);
    
end

legend(ax,unq);


end

