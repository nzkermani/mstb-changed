function [ topSc,mz,finalSpec ] = rawREIMSfolder( folder )
%rawREIMSfolder - read in a bunch of files from a single folder

if nargin == 0
    folder = 'E:\Dropbox\Imperial\Projects\OMB\for_james\naocl_masslynx_processed\';
end

% Find the files...
files = fileFinderAll(folder,'raw',true);
numF = size(files,1);

% Somewhere to save the most intense scans
topSc = cell(numF,1);
allSc = cell(numF,1);

% Now just look through...
for n = 1:numF
    
    % Run the function...
    %[mz,x,xy,opts] = rawREIMS([files{n,1} filesep files{n,2}],opts);
    [sp,tot] = rawREIMSread([files{n,1} filesep files{n,2}]);
    
    % Bin the data...
    [mz,sp1] = dbBinningFixed(sp(1),[50 1200],0.01);
    allDat = zeros(size(sp,2),numel(mz));
    for r = 1:size(sp,2)
        [~,tmp] = dbBinningFixed(sp(r),[50 1200],0.01);
        allDat(r,:) = tmp.al;
    end
    
    % If n == 1
    if n == 1
        finalSpec = zeros(numF,numel(mz));
    end
    finalSpec(n,:) = nanmean(allDat,1);
        
    % Find the most intense scan
    [~,scanIdx] = max(tot);
    
    % Keep this one...
    %topSc{n,1} = sp{scanIdx};
    
    % Keep all the data?
    %allSc{n,1} = sp;
        
end

end

