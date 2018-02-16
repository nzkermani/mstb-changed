function [mz,data] = dbBinningFixed(data,mzRange,res)
%dbBinningFixed - perform binning of the data according to the resolution
%provided and over the m/z range specified

numF = size(data,2);

% Create the m/z vector according to the resolution
mzMin = min(mzRange);
mzMax = max(mzRange);
mz = mzMin:res:mzMax;

% Now that we have the resolution, we can start the process
for n = 1:numF
    
    % Trim data according to the m/z range
    tmpMZ = data(n).mz;
    mask = tmpMZ >= mzMin & tmpMZ <= mzMax;
    tmpMZ = tmpMZ(mask);
    tmpSP = data(n).sp(:,mask);
    
    % From the file's m/z vector, determine the indices
    idx = floor(tmpMZ / res) - (mzMin / res) + 1;
        
    % Now create an empty matrix of the new size
    data(n).al = zeros(size(tmpSP,1),numel(mz));
    
    % Add in the old variables
    for r = 1:numel(idx)
        data(n).al(:,idx(r)) = data(n).al(:,idx(r)) + tmpSP(:,r);
    end
    
    if sum(data(n).al(:)) ~= sum(tmpSP(:))
        %disp('Summation error');        
        try
            disp(data(n).file);
        catch
        end            
    end

end

mz = mz';
    
end

