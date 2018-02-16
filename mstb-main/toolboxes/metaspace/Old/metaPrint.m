function [ output_args ] = metaPrint( heads,meta )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

numF = numel(heads);


fid = fopen('/Users/jmckenzi/Dropbox/MTSP-Meta.csv','w');

for n = 1:numF
    
    fprintf(fid,'%s,',heads{n});
    
    % Determine unique entries
    unq = unique(meta(:,n));
    numU = numel(unq);
    
    for r = 1:numU
        
        fprintf(fid,'%s,',unq{r});
        
    end
    
    fprintf(fid,'\n');
    
end

end

