function [ op,days ] = day2num(tt)
%day2num - convert Monday, Tuesday, etc to a series from 1-7

days = {'Monday';'Tuesday';'Wednesday';'Thursday';'Friday';'Saturday';'Sunday'};

op = zeros(size(tt));

for n = 1:7
    
    fx = strcmpi(tt,days{n});
    
    op(fx) = n;
    
end

end

