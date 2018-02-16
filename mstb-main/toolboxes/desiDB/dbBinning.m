function [mz,data] = dbBinning(data,res)
%dbBinning - perform binning of the data according to the resolution
%provided

numF = size(data,2);

% Determine the min and max m/z values across all samples
vals = zeros(numF,2);
for n = 1:numF
    vals(n,:) = [min(data(n).mz) max(data(n).mz)];
end

% Create the m/z vector according to the resolution
mzMin = floor(min(vals(:,1)));
mzMax = ceil(max(vals(:,2)));
mz = mzMin:res:mzMax;

% Now that we have the resolution, we can start the process
for n = 1:numF
    
    % From the file's m/z vector, determine the indices
    tmpMZ = data(n).mz;    
    idx = floor(tmpMZ / res) - (mzMin / res) + 1;
    
    % Now create an empty matrix of the new size
    data(n).al = zeros(size(data(n).sp,1),numel(mz));
    
    % Add in the old variables
    for r = 1:numel(idx)
        data(n).al(:,idx(r)) = data(n).al(:,idx(r)) + data(n).sp(:,r);
    end
    
    if sum(data(n).al(:)) ~= sum(data(n).sp(:))
        disp('Summation error');        
        try
            disp(data(n).file);
        catch
        end            
    end

end

mz = mz';
    
end

