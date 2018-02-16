function [allv] = detVar(resSt,vec)
%detVar - determine the variances of the data, but it requires
%concatenation beforehand
%
% resSt is an (n,m) structure where:
% n = number of scans
% m = number of interpolation methods
%
% So we concatenate over resSt(:,m) for comparison of variances

[valN,valM] = size(resSt)

% Structure for storing the variances
allv = struct('mz',[],'var',[],'cs',[]);

for m = 1:valM
    
    tmp = vertcat(resSt(:,m).sp);
    
    %allv(m).mz  = vec(m).mz;
    allv(m).var = var(tmp,[],1);
    
    allv(m).cs  = cumsum(allv(m).var);
    
end

% figure; hold on;
% h = zeros(valM,1);
% cols = jet(valM);
% for m = 1:valM
%     
%     h(m,1) = plot(allv(m).mz,100 .* allv(m).cs ./ max(allv(m).cs),...
%         'Color',cols(m,:),...
%         'LineWidth',2);
%     
% end



end

