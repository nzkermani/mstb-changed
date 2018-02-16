function [ sm ] = movingMZWindow(mz,sp,da,method)
%movingWindow - apply the function in 'method' across the spectrum
%considering points from -hw:0:+hw
%
% James McKenzie, 2016

numV = numel(sp);
sm = zeros(size(sp));

switch method
    case 'min'
        smFun = @(x) min(x);
    case 'mean'
        smFun = @(x) mean(x);
    case 'median'
        smFun = @(x) median(x);
    case 'max'
        smFun = @(x) max(x);
    case '10'
        smFun = @(x) prctile(x,10);
    otherwise
        smFun = @(x) x;
end

for n = 1:numV
    
    % Define the boundaries over which to smooth
    mzl = max([min(mz) mz(n)-da]);
    mzh = min([mz(n)+da max(mz)]);
    
    fx = mz >= mzl & mz <= mzh;
    
    sm(n) = smFun(sp(fx));

% %     % Smooth...
% %     switch method
% %         
% %         case 'min'
% %             sm(n) = min(sp(fx));
% %             
% %         case 'mean'
% %             sm(n) = mean(sp(fx));
% %             
% %         case 'median'
% %             sm(n) = median(sp(fx));
% %             
% %         case 'max'
% %             sm(n) = max(sp(fx));
% %             
% %         case '10'
% %             sm(n) = prctile(sp(fx),10);
% %             
% %         otherwise
% %             sm(n) = sp(n);
% %     end
end

end

