function [ pvar,pq ] = i4iSS(data)
%i4iSS - determine the sample size for the i4i data

warning off all

% Determine the mean etc values for the spectra
[m1,m2,s1,s2,ss] = calcMeanSTD(data.sp,data.meta.class);

% For each variable, determine its power...
numV = size(data.sp,2);
pvar = zeros(numV,2);
pq = ones(numV,2);
[unq,~,ind] = unique(data.meta.class);

ratio = sum(ind == 1) / sum(ind == 2)
if ratio < 1
    ratio = 1 / ratio;
end

for n = 1:numV
%     pvar(n,1) = sampsizepwr('t2',...
%         [m1(n) s1(n)],m2(n),...
%         [],min([sum(ind == 1) sum(ind == 2)]),...
%         'Alpha',0.01,...
%         'Ratio',ratio,...
%         'Tail','both');
    
    pvar(n,2) = sampsizepwr('t2',...
        [m1(n) s1(n)],m2(n),...
        0.95,[],...
        'Alpha',0.01,...
        'Ratio',ratio,...
        'Tail','both');
    
    % Calculate the p value using KW
    pq(n,1) = anova1(data.sp(:,n),data.meta.class,'off');
    
    if mod(n,100) == 0
        disp(int2str(n));
    end
end

% Convert the p values to q values
pq(:,2) = p2q(pq(:,1),0.05);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m1,m2,s1,s2,ss] = calcMeanSTD(sp,grp)

[~,~,ind] = unique(grp);
i1 = ind == 1;
i2 = ind == 2;

m1 = nanmean(sp(i1,:),1);
m2 = nanmean(sp(i2,:),1);
s1 = std(sp(i1,:),[],1);
s2 = std(sp(i2,:),[],1);
ss = std(sp,[],1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%