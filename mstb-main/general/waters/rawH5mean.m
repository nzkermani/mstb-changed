function [ mz,sp,pmz ] = rawH5mean(file,mzRange)
%rawH5mean - determine mean spectrum

% Define mz values if not provided
if nargin == 1
    mzRange = [100 1000];
end

% Determine number of pixels
totInt = h5read(file,'/totPts');
numS = numel(totInt);
    
% Paolo's code
pmz = cell(numS,1);

wb = waitbar(0,'Hello');

% Loop through each pixel
for n = 1:numS
    
    % Scan name
    sn = ['/raw/scan/' int2str(n)];
    
    % Read spectrum...
    tmp = h5read(file,sn);
    mask = tmp(:,1) >= min(mzRange) & tmp(:,1) <= max(mzRange);
    
    % Paolo's code
    pmz{n} = tmp;

    
    if n == 1
        all = tmp(mask,:);
    else
        
        % Concat the two spectra, sort, and then combine the duplicated mz
        % values
        all = cat(1,all,tmp(mask,:));
        
        all = sortrows(all,1);
        
        % Determine unique values...
        [a,~,c] = unique(all(:,1));
        
        % Accumulate the duplicated values.
        aa = accumarray(c,all(:,2));
        
        % Save...
        all = [a aa];       
        
        
    end
    
    waitbar(n/numS,wb);
end

mz = all(:,1);
sp = all(:,2) / numS;

% Paolo's code
pmz = cell2mat(pmz);
[unqMZ,~,unqInd] = unique(pmz);

delete(wb);

end

