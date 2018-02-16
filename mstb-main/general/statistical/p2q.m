function [q] = p2q(p,alpha)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
    alpha = 0.05;
end

% Check that is a vector
sz = size(p);
if min(sz) ~= 1
    error('Need a vector of p values');
end

% Transpose if the wrong way round
if sz(1) > sz(2)
    p = p';
end

% Use this method
[q,~] = getBHYqVls(p,alpha);

end

