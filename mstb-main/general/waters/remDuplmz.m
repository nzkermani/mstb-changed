function [mz, counts] = remDuplmz(mz, counts)
% Remove duplicated m/z values after rounding of according the tolerance of an instrument   
diffmz            = diff(mz)';
[mzuniq, mzindcs] = find(diffmz>1);
dubmzindcs        = find(diff(mzindcs)>1);
if ~isempty(dubmzindcs)
    counts_temp = counts([1 mzindcs+1]);
    for i = dubmzindcs
        counts_temp(i+1) = sum(counts(mzindcs(i)+1:mzindcs(i+1)));
    end
    mz     = mz([1 mzindcs+1]);
    counts = counts_temp;
end