function [ xp,yp ] = insertZeros( x,y,s,plotFlag )
%insertZeros - run prior to plotting, inserts zeros into spectra between
%variables, to prevent stupidly broad variables exemplified by plot(x,y);
% x         - the m/z values etc...
% y         - spectral variables
% s         - the spacing, which should ideally be tiny
% plotFlag  - true/false for plotting the spectra

if nargin == 2
    s = [0.01 0.01];
    plotFlag = false;
elseif nargin == 3
    % Decide if this was a logical thing or a number?
    ws = whos('s');
    if strcmp(ws.class,'logical');
        plotFlag = true;
        s = [0.01 0.01];
    else
        plotFlag = false;
        s = [0.01 0.01];
    end
end

% Have it for each side...
if numel(s) == 1
    s = [s s];
end

% How many spectral variables are there?
numV = numel(x);

% Create a new vector that is thrice the size of the original vector...
xp = zeros(1,numV*3);
yp = zeros(size(y,1),numV*3);

% Create the map, i.e. 1_2_3 4_5_6 7_8_9 where 2,5,8 are the peaks and
% others are the zero intensities at ±s
map = 2:3:(numV*3)-1;

% Transfer the mapping
xp(1,map) = x;
yp(:,map) = y;

% Now need to add in mz values for each at either side of the main peak
mapL = 1:3:(numV*3)-2;
xp(1,mapL) = x - s(1);

mapR = 3:3:(numV*3);
xp(1,mapR) = x + s(2);%s;

if plotFlag
    figure;
    plot(xp,yp);
end

return

% Delete overlapping variables...
ddd = diff(xp);
ddd = find(ddd < 0);
ddd = ddd + 1;

xp(ddd) = [];
yp(ddd) = [];

end

