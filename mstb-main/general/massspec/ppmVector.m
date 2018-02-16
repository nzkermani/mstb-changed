function [ vec ] = ppmVector(mzL,mzH,ppm)
%ppmVector - create an m/z vector with ppm dependent gaps. The result is
%non-linear spacing between the m/z values
%
% This is V2, which now accepts polynomial functions to determine the ppm
% spacing. This is a marked improvement, see interpExp for comparision...

% Can we estimate the size of the matrix required? This will always produce
% a much larger vector than necessary
if numel(ppm) == 1
    i = (mzH-mzL) ./ (ppm * [mzL mzH] / 1e6);
else
    % Then we have provided a fitted function for the ppm values to
    % increase rather than by a constant amount each time. Thus this gets a
    % lot more complicated!
    meanMZ = (mzH-mzL)/2;    
    meanPPM = polyval(ppm,meanMZ);    
    i = (mzH-mzL) ./ (meanPPM * [mzL mzH] / 1e6);
end

% Common to that above...
maxP = round(mean(i));
vec = zeros(maxP,1);    
        
% Let's go
n = 1;
vec(n,1) = mzL;

% Loop
while vec(n,1) < mzH
    
    % Calculate x ppm of this mz value
    if numel(ppm) == 1
        i = ppm * vec(n,1) / 1e6;
    else
        i = polyval(ppm,vec(n,1)) * vec(n,1) / 1e6;
    end

    % Increment counter...
    n = n + 1;
    
    % ...add in the value
    vec(n,1) = vec(n-1,1) + i;
    
end

% Trim...
vec = vec(1:n,1);

end

