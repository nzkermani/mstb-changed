function [orig] = histoLoad
%histoLoad - load a histological image

% Default path
imgPath = ['/Users/jmckenzi/DB/'];

% Get a file
[imgName,imgPath] = uigetfile({'*.png; *.jpeg; *.jpg; *.tiff'},...
    'Select',imgPath);

orig = imread([imgPath imgName]);

end

