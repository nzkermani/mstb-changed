function [ mask ] = statsVariableFiltering(sts)
%statsVariableFiltering - perform filtering of variables according to
%certain conditions...

% Need to determine how many files there are, as this is to be done on a
% file by file basis, rather than with reference to observations...
[unqF,firstF,indF] = unique(sts.raw.meta.fileID);
numF = numel(unqF);

% Store variable presence
pre = zeros(numF,size(sts.raw.sp,2));

% And dates
dates = zeros(numF,1);

% Now run through obs from each file and
for n = 1:numF
    
    fx = indF == n;
    
    pre(n,:) = sum(sts.raw.sp(fx,:),1) > 0;
    
    dates(n,1) = sts.raw.meta.date(firstF(n));
    
end

% How about filtering according to appearing in half of the files
frq = nansum(pre > 0,1);

mask = frq > (size(pre,1) * 0.5);

return


% Histogram the number of variables according to acquisition date
sumV = sum(pre,2);

% For each observation, write out the number of variables per file
totV = zeros(size(indF));
for n = 1:numF
    
    fx = indF == n;
    totV(fx) = sumV(n,1);
end

sts.raw.meta.totV = totV;
sts.raw.meta.numV = nansum(sts.raw.sp > 0,2);

    
end

