function [extr,rawmz,totsum] = rawImageMake(sp,xy2D,mz,ppm)
%rawImageMake - need to make capable to handle multiple ion inputs.
%
% LMC has been removed

% Determine ppm tolerances...
tol = ppm * mz / 1e6;

% How many ions are there to be extracted?
numI = numel(mz);

% Size of the image?
sz = size(xy2D);

% Storage
extr = NaN(sz(1),sz(2),numI);
rawmz = NaN(size(extr));
totsum = NaN(sz(1),sz(2));

% Loop
for x = 1:sz(1)
    for y = 1:sz(2)
    
        % Image coordinate
        n = xy2D(x,y);
        if isnan(n)
            continue;
        end
        if isempty(sp{n})
            continue;
        end
        
        % Loop through for each ion
        for i = 1:numI
    
            % Find the biggest peak within ppm ppm
            mask = sp{n}(1,:) >= (mz(i)-tol(i)) & sp{n}(1,:) <= (mz(i)+tol(i));
            tmp = sp{n}(2,:) .* mask;
            [val,idx] = max(tmp);
            if val > 0
                extr(x,y,i) = sp{n}(2,idx);
                rawmz(x,y,i) = sp{n}(1,idx);
            end            
        end    
        
        % Total sum
        totsum(x,y) = nansum(sp{n}(2,:));
        
    end    
end

end

