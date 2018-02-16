function [ output_args ] = multiAnova( mz,sp,meta,names )
%multiAnova - can we identify variables which fail to be useful?

numV = size(sp,2);
numG = size(meta,2);
pval = zeros(numV,numG);

for n = 1:numV
    
    [p,t,s] = anovan(sp(:,n),meta,...
        'Display','off',...
        'VarNames',names);
    
    pval(n,:) = p';
    
end



end

