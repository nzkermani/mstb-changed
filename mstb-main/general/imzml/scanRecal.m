function [ idx ] = scanRecal(mz,sp,ions,hw)
%scanRecal - recalibrate a single scan using some diagnostic ions...

% Use these ions as the fixed points
% ions = [255.2329 281.2485 465.3043 536.5047 744.5548 885.5498];
% hw = 20;

% Find the ions in the scan...
idx = NaN(size(ions));

for n = 1:numel(ions)
    
    % Find largest signal - this will be problematic when we talk about the
    % background...
    df = hw * ions(n) / 1e6;
    fx = mz > (ions(n) - df) & mz < (ions(n) + df);
    if sum(fx) > 0
        fx = fx ./ sp;    
        [~,idx(n)] = max(fx);
    end
end

idx(~isnan(idx)) = mz(idx(~isnan(idx)));

return


fx = ~isnan(idx);


% Determine the differences between the ions
df = bsxfun(@minus,ions(fx),mz(idx(fx))');

p = polyfit(ions(fx),df,3);

y1 = polyval(p,mz);


figure; hold on;
scatter(ions(fx),mz(idx(fx)));
plot(ions,y1);

end

