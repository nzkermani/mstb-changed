function [ output_args ] = mtspCompareOldNew(old,new,recal)
%mtspCompareOldNew - compare old and new annoations quantities

% There may be disparities in the number of files of each one, so just try
% to focus on the new files, where present...
numN = size(new.group,1);

% Store the old|new|ppm deviation in here
match = zeros(numN,3);

% Loop through
for n = 1:numN
    
    % This file name
    tf = new.group{n,1}(1:end-10)
    
    % Match to the old
    fx = strcmp(old.group(:,1),[tf '.mat']);
        
    % Combine annotations
    if sum(fx) == 1
        match(n,1:2) = [old.numAnno(fx,1) new.numAnno(fx,1)];
    end
    
    
end




end

