function [op] = fitGauss(xdata,y)
% Fit the gauss function to the data provided. Method taken from the
% website: http://stackoverflow.com/questions/13290508/...
% how-to-fit-a-gaussian-to-data-in-matlab-octave (w/o the dots)
%figure; stem(xdata,y);

% Design matrix for least squares fit
xdata = xdata(:);
A = [xdata.^2 xdata ones(size(xdata))];

% log(y)
b = log(y(:));
%cb = isinf(b);
%b(cb) = 0;

% Least squares solution for x
x = A \ b;

% Calculate gaussian parameters
mu = -x(2)/x(1)/2;
sigma = sqrt(-1/2/x(1));

% What is the indefinite integral of the function?
maxI = max(y);
ii = maxI * sigma * sqrt(2 * pi);

% Here is the output: mu/sigma/intensity/integral
op = [mu sigma maxI ii];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
