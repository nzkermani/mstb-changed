function [ output_args ] = mtspDetermineNumAnnos
%mtspDetermineNumAnnos - read in the files to determine annotation
%quantities

path = '/Users/jmckenzi/Dropbox/Imperial/Projects/Metaspace/Engine Dump/MTSP';
allF = fileFinderAll(path,'mat',true);
numF = size(allF,1);

sz = zeros(numF,1);

for n = 1:numF
    
    load([allF{n,1} filesep allF{n,2}],'mz');
    
    sz(n,1) = numel(mz);
    
end
        


end

