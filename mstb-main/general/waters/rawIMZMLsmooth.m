function [ output_args ] = rawIMZMLsmooth(mz,avg,frq)
%rawIMZMLsmooth - simple smoothing to get peak shapes rather than jagged
%things

sma = zeros(size(avg));
for n = 1:size(avg,2)
    
    sma(:,n) = sgolayfilt(frq(:,n),5,21);
    
end


end

