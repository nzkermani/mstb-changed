function [ sm ] = movingWindow(sp,hw,method)
%movingWindow - apply the function in 'method' across the spectrum
%considering points from -hw:0:+hw
%
% James McKenzie, 2016

numV = numel(sp);
sm = zeros(size(sp));

% Determine the function out here
switch method

    case 'min'
        myfun = @min;%(sp(st:fn));

    case 'mean'
        myfun = @mean;%(sp(st:fn));

    case 'median'
        myfun = @median;%(sp(st:fn));

    case 'max'
        myfun = @max;%(sp(st:fn));

    case '10'
        myfun = @prctile;%(sp(st:fn),10);
        
    case 'sum'
        myfun = @sum;
        
    case 'gauss'
        myfun = @jsmGauss;

    otherwise
        error('pointless');        
end


for n = 1:numV
    
    % Define the boundaries over which to smooth
    st = max([1 n-hw]);
    fn = min([n+hw numV]);
    
    % Smooth...
    sm(n) = myfun(sp(st:fn));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sm] = jsmGauss(sp)
% Create a gauss filter

xe


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%