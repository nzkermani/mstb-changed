function [res] = mtspCSVmultiple
%mtspCSVmultiple - etc etc.

% Define the path
path = '/Users/jmckenzi/Dropbox/Imperial/Projects/Metaspace/ICL/';

% Find all CSV files in the folder
csvs = fileFinderAll(path,'csv');
csvs(1,:) = [];
numF = size(csvs,1);
res = cell(numF,1);

% Waitbar
wb = waitbar(0);

% Loop through
for n = 1: numF
    
    [~,res{n,1}] = mtspCSVimport([csvs{n,1} filesep csvs{n,2}]);
    
    waitbar(n/numF,wb,[int2str(n) '/' int2str(numF)]);
end

delete(wb);


end

