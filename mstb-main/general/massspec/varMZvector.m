function [ vec ] = varMZvector(mzL,mzH,fit)
%varMZvector - create an m/z vector with ppm dependent gaps. The result is
%non-linear spacing between the m/z values
%
% This is V3, which only accepts polynomials and for non-PPM differences as
% these are non-linear
if numel(fit) == 1
    error('Wrong size polynomial fit');
end

% Can we estimate the size of the matrix required? This will always produce
% a much larger vector than necessary
meanMZ = (mzH-mzL)/2;    
meanDiff = mean([polyval(fit,mzL) polyval(fit,mzH)]);
i = (mzH-mzL) ./ meanDiff;

% Common to that above...
maxP = round(mean(i));
vec = zeros(maxP,1);    
        
% Let's go
n = 1;
vec(n,1) = mzL;

% Loop
while vec(n,1) < mzH
    
    % Calculate increment of this value
    i = polyval(fit,vec(n,1));
    
    % Increment counter...
    n = n + 1;
    
    % ...add in the value
    vec(n,1) = vec(n-1,1) + i;
    
end

% Trim...
vec = vec(1:n,1);

end

