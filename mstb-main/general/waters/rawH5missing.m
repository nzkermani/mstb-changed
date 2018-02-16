function [ output_args ] = rawH5missing(allF)
%rawH5missing - find missing files, which need to be transferred or
%processed

% Where the final files are kept
procPath = '/Users/jmckenzi/Documents/Box Sync/H5 Convert/';

h5Files = fileFinderAll(procPath,'h5');

numF = size(allF,1);
chk = false(numF,1);

for n = 1:numF
    
    % Which file to find?
    ffind = allF{n,1};
    
    % Replace .raw with .h5
    dot = strfind(ffind,'.');
    ffind = [ffind(1:dot(end)) 'h5'];
        
    % Match to h5Files
    cmp = strcmp(h5Files(:,2),ffind);
    
    if sum(cmp) == 1
        chk(n,1) = true;
    end
    
end




end

