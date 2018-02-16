function [ output_args ] = interpShift(mz,sp,interpMZ,res)
%interpShift - shift a spectrum by a specified amount to see how the
%intepolated data varies...

shift = [-1:0.1:1];
numS = numel(shift);

spec = zeros(numS,numel(interpMZ));

for n = 1:numS
    
    % Modify the mz vector by shift(n) ppm
    ppm = shift(n) * mz / 1e6;    
    newMZ = mz + ppm;
    
    % Do the interpolation
    spec(n,:) = interp1(newMZ,sp,interpMZ,'linear');
    
end

figure; plot(interpMZ,spec);
title(res);

figure; scatter(interpMZ,std(spec,[],1)./mean(spec,1));
title(res);

end

