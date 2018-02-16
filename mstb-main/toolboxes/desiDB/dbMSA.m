function [data,files] = dbMSA(files,opts)
%dbMSA - DESI multiple sample analysis, which is launched through the
%desiDB function. We need to read in the files, align the m/z vectors, and
%then save the results. Everything else happens elsewhere, by loading the
%saved results and then working some magic.

% Import the files - this is different for Metaspace and normal approaches
switch opts.pixel
    case 'Annotated'
        [data,files] = importFiles(files,opts);
    case 'Metaspace'
        
        % Read in the data
        [data] = importMetaspace(files,opts);
        
        % Need to align the annotations...
        [data] = mtspAlignAnnotations(data);
        assignin('base','data',data);
        assignin('base','opts',opts);
        %files = [];
        return
end

assignin('base','data',data);
assignin('base','opts',opts);

% Run the alignment
if opts.mzRes < 1 && ~isnan(opts.ppmTol)
    %[data,cmz] = mzAlignment(data,opts);
    [data,cmz,qcmz] = mzAlign(data,opts.ppmTol);

elseif isnan(opts.ppmTol)
    % Perform binning - but don't tell anyone!
    [cmz,data] = dbBinning(data,opts.mzRes);    
end

% Combine the spectra together, and make vectors so that we know which file
% each spectrum belongs to.
numF = size(data,2);
allInfo = struct('date',[],'fileID',[],'histID',[],'foldID',[],'foldID2',[],'foldID3',[]);
for n = 1:numF
    tmp = data(n).file;
    dot = strfind(tmp,'.');
    
    allInfo(n).date   = repmat(data(n).date,[size(data(n).al,1) 1]);
    allInfo(n).fileID = repmat({tmp(1:dot(end)-1)},[size(data(n).al,1) 1]);
    allInfo(n).histID = data(n).histID;
    allInfo(n).foldID = repmat({data(n).subType},[size(data(n).al,1) 1]);
    
    try
        [aa,bb] = previousFolder(data(n).path);
        [cc,dd] = previousFolder(aa);
        [ee,ff] = previousFolder(cc);
        allInfo(n).foldID2 = repmat({dd},[size(data(n).al,1) 1]);
        allInfo(n).foldID3 = repmat({ff},[size(data(n).al,1) 1]);
    catch
        allInfo(n).foldID2 = repmat({'Fail'},[size(data(n).al,1) 1]);
        allInfo(n).foldID3 = repmat({'Fail'},[size(data(n).al,1) 1]);
    end
    
end
full.mz = cmz;
full.sp = vertcat(data.al);
full.meta.date   = vertcat(allInfo.date);
full.meta.histID = vertcat(allInfo.histID);
full.meta.fileID = vertcat(allInfo.fileID);
full.meta.foldID = vertcat(allInfo.foldID);
full.meta.foldID2 = vertcat(allInfo.foldID2);
full.meta.foldID3 = vertcat(allInfo.foldID3);
data = full;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,files] = importFiles(hits,opts)
% Read in the pertinent parts of the files

% Waitbar
wb = waitbar(0,'Import files');

% Data storage
numF = size(hits,1);
data = struct('path',[],'file',[],'date',[],'raw',[],'subType',[],...
    'mz',[],'sp',[],'histID',[],'al',[]);
skip = true(numF,1);

% Loop through the files
for n = 1:numF
    
    % Read in each one
    tmp = open([hits{n,1} filesep hits{n,2}]);
    data(n).path = hits{n,1};
    data(n).file = hits{n,2};
    
    % Check that this file is the right type
    if ~isfield(tmp,'dpn');
        skip(n,1) = false;
        continue;
    end
    
    if ~isfield(tmp.dpn,'anno');
        skip(n,1) = false;
        continue
    end
    
    % Also check that files are within the date period set. But if we have
    % no date in the file then we have to include it...
    if isfield(tmp.dpn,'date')
        
        if isempty(tmp.dpn.date)
            tmp.dpn.date = mean([opts.dateLow opts.dateHigh]);
        end
        
        % Check that it lies between the values...
        if tmp.dpn.date >= opts.dateLow & tmp.dpn.date <= opts.dateHigh
            % Keep the file
        else
            skip(n,1) = false;
            continue;
        end
        
    end
    
    % Save the raw file name...
    data(n).raw = [tmp.dpn.file.dir tmp.dpn.file.nam];
    
    % If there is a date field, then take it...
    if isfield(tmp.dpn,'date')
        data(n).date = tmp.dpn.date;
    else
        data(n).date = now;
    end
    
    % Determine the folder (i.e. tissue type) of the file... And then
    % create an array of the same same as histID
    sl = strfind(data(n).path,filesep);
    data(n).subType = data(n).path(sl(end)+1:end);

    % Perform the annotation consolidation function!
    [tmp.dpn.anno] = xxxAnnotationConsolidate(tmp.dpn.anno);
    
    % Determine the histologically annotated pixels
    [isAnno,histID,~] = desiAnnotationExtract(tmp.dpn);
    
    
    % Now we need to determine how many pixels are in each tissue type? We
    % do this for each tissue type individually
    [unqA,~,~] = unique(isAnno);
    [frqA,~] = hist(isAnno,unqA);
    
    % Find any frqA with values less than opts.minPix
    for f = 1:numel(frqA)
        if frqA(f) < opts.minPix
            fx = isAnno == unqA(f);
            isAnno(fx) = 0;
            %histID(fx) = '';
        end
    end
    
    % Which pixels remain?
    mask = isAnno ~= 0;
    
    % Skip if there are no pixels...
    if sum(isAnno(mask)) == 0
        skip(n,1) = false;
        continue;
    end
    
    % Save the MS information
    mzVec = tmp.dpn.d1.mz;
    
    % Trim according to m/z range desired
    mzMask = mzVec >= opts.mzRange(1) & mzVec <= opts.mzRange(2);    
    
    % Reshape the matrix
    sz = size(tmp.dpn.d1.sp);
    sp = reshape(tmp.dpn.d1.sp,[sz(1)*sz(2) sz(3)]);
    
    % Trim bits out of the mz and sp vector/matrix
    data(n).mz = mzVec(mzMask);
    if size(data(n).mz,1) > size(data(n).mz,2)
        data(n).mz = data(n).mz';
    end
    sp = sp(:,mzMask);
    
    % If this comes from old data, then we need to consider that it might
    % be logged
    if max(sp(:)) < 100 && min(sp(:)) > 0
        mn = min(sp(:));        
        sp = exp(sp) - exp(mn);        
    end
    
    % Warning in case of issues
    if min(sp(:)) ~= 0
        disp('Potential issues with non-zero baseline');
    end
    
    % Here is the place to decide which of the annotated pixels we want.
    % Where 'all' is selected, then we need to do exactly nothing.
    % Otherwise, we have to average or subsample.
    switch opts.pixelCombo        
        case 'All'
            data(n).sp = sp(mask,:);
            data(n).histID = histID(mask,:);
            
        case 'All-QC'
            
            % Determine unique tissue types in this file...
            [unqTT,~,indTT] = unique(isAnno);
            newMask = false(size(mask));
            for p = 2:numel(unqTT)
                
                % These pixels
                px = indTT == p;
                ax = histID{find(px == 1,1,'first')};
                
                % Is this the background?
                isBG = any(strfind(lower(ax),'background')) || any(strfind(lower(ax),'bg'));
                if isBG
                    [inc] = specCheck(sp(indTT == p,:),true);
                else
                    [inc] = specCheck(sp(indTT == p,:));
                end
                
                % How to reintegrate into 'mask' which determines which
                % pixels we include
                newMask(px) = inc;
            end
            
            % Now we need to save the observations that passed the tests...
            
            % Perhaps we can PCA of all pixels and just the selected ones
            % to see if there is any difference...
            [~,sa,~] = pca(sp(mask,:));
            [~,sb,~] = pca(sp(newMask,:));
            
            figure;
            ax1 = subplot(1,2,1);
            [ca,~,ia] = unique(histID(mask,:));
            scatter(sa(:,1),sa(:,2),80,ia,'o','filled');
            
            ax2 = subplot(1,2,2);
            [cb,~,ib] = unique(histID(newMask,:));
            scatter(sb(:,1),sb(:,2),80,ib,'o','filled');
            
            
            % Save...
            data(n).sp = sp(newMask,:);
            data(n).histID = histID(newMask,:);            
            
        otherwise
            % This requires some portion of the pixels
            sp = sp(mask,:);
            histID = histID(mask);
            [grp,~,ind] = unique(histID);
            numG = numel(grp);
            
            % Create somewhere to store the data...
            switch opts.pixelCombo
                case {'TT Mean','TT Median'}
                    tmpSp = zeros(numG,size(sp,2));
                case 'Random 5'
                    pickMe = zeros(5,numG);
            end
            
            % Loop through each tissue type
            for r = 1:numG
                
                % Indices of these pixels
                fx = ind == r;
                
                % Now we need to decide how to treat them
                switch opts.pixelCombo
                    case 'TT Mean'
                        tmpSp(r,:) = nanmean(sp(fx,:),1);
                    case 'TT Median'
                        tmpSp(r,:) = nanmedian(sp(fx,:),1);
                    case 'Random 5'                        
                        fy = find(fx);
                        numRand = min([5 numel(fy)]);
                        fz = randperm(numel(fy),numRand);
                        pickMe(1:numRand,r) = fy(fz);
                end
            end
            
            % Once we've looped through then we need to finish off
            switch opts.pixelCombo
                case 'Random 5'
                    pickMe = pickMe(:);
                    pickMe = pickMe(pickMe ~= 0);
                    data(n).sp = sp(pickMe,:);
                    data(n).histID = histID(pickMe,:);
                    
                case {'TT Mean','TT Median'}
                    data(n).sp = tmpSp;
                    data(n).histID = grp;
            end
            

    end
        
    % Update
    waitbar(n/numF,wb);
end

delete(wb);

% Trim out fake files
data = data(skip);
files = hits(skip,:);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = importMetaspace(hits,opts)
% Read in the pertinent parts of the files

% Waitbar
wb = waitbar(0,'Import files');

% Data storage
numF = size(hits,1);
data = struct('path',[],'file',[],'subType',[],'annos',[],'sp',[],'histID',[],'al',[]);

% Loop through the files
for n = 1:numF
    
    % Read in each one
    tmp = open([hits{n,1} filesep hits{n,2}]);
    data(n).path = hits{n,1};
    data(n).file = hits{n,2};
    
    % Check that this file actually contains MTSP information
    if isfield(tmp,'dpn')        
        if ~isfield(tmp.dpn,'mtsp');
            disp(['NO MTSP: ' data(n).file]);
            continue;
        end
    else
        continue;
    end
    
    % Determine the folder (i.e. tissue type) of the file... And then
    % create an array of the same same as histID
    sl = strfind(data(n).path,filesep);
    data(n).subType = data(n).path(sl(end)+1:end);

    % Perform the annotation consolidation function!
    [tmp.dpn.anno] = xxxAnnotationConsolidate(tmp.dpn.anno);
    
    % Determine the histologically annotated pixels
    [isAnno,histID,~] = desiAnnotationExtract(tmp.dpn);
    mask = isAnno ~= 0;

    % Need to check the sizes of the matrices - they should be identical,
    % but perhaps consider shaving off the bottom row...
    szOrig = size(tmp.dpn.d1.sp);
    szMtsp = size(tmp.dpn.mtsp.sp);
    
    % If the number of rows don't agree, trim off the bottom...
    sp = tmp.dpn.mtsp.sp(1:szOrig(1),1:szOrig(2),:);
    sz = size(sp);
    if numel(sz) == 2
        sz = [sz 1];
    end
    
    % Reshape the MTSP matrix
    sp = reshape(sp,[sz(1)*sz(2) sz(3)]);
    
    % Add in the annotations - note these need to replace the m/z values
    data(n).annos = tmp.dpn.mtsp.annos;
    
    % Here is the place to decide which of the annotated pixels we want.
    % Where 'all' is selected, then we need to do exactly nothing.
    % Otherwise, we have to average or subsample.
    switch opts.pixelCombo        
        case 'All'
            data(n).sp = sp(mask,:);
            data(n).histID = histID(mask,:);
            
        case 'All-QC'
            error('Not compatible here');
                        
        otherwise
            % This requires some portion of the pixels
            sp = sp(mask,:);
            histID = histID(mask);
            [grp,~,ind] = unique(histID);
            numG = numel(grp);
            
            % Create somewhere to store the data...
            switch opts.pixelCombo
                case {'TT Mean','TT Median'}
                    tmpSp = zeros(numG,size(sp,2));
                case 'Random 5'
                    pickMe = zeros(5,numG);
            end
            
            % Loop through each tissue type
            for r = 1:numG
                
                % Indices of these pixels
                fx = ind == r;
                
                % Now we need to decide how to treat them
                switch opts.pixelCombo
                    case 'TT Mean'
                        tmpSp(r,:) = nanmean(sp(fx,:),1);
                    case 'TT Median'
                        tmpSp(r,:) = nanmedian(sp(fx,:),1);
                    case 'Random 5'                        
                        fy = find(fx);
                        numRand = min([5 numel(fy)]);
                        fz = randperm(numel(fy),numRand);
                        pickMe(1:numRand,r) = fy(fz);
                end
            end
            
            % Once we've looped through then we need to finish off
            switch opts.pixelCombo
                case 'Random 5'
                    pickMe = pickMe(:);
                    pickMe = pickMe(pickMe ~= 0);
                    data(n).sp = sp(pickMe,:);
                    data(n).histID = histID(pickMe,:);
                    
                case {'TT Mean','TT Median'}
                    data(n).sp = tmpSp;
                    data(n).histID = grp;
            end
            

    end
    
        
    % Update
    waitbar(n/numF,wb);
end

delete(wb);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,cmz] = mzAlignment(data,opts)
% First define the cmz vector, and then align everything to it

% Define the peak matching options
opts.estimationMethod = 'histogram';
opts.display = false;
opts.handBand = eval(['@(z) (' num2str(opts.ppmTol) '* z / 1e6)']);            
opts.handGap  = eval(['@(z) (' num2str(opts.ppmTol) '* z(:,1) / 1e6)']);
opts.handDist = @(R,S) abs(sum((R-S),2));
opts.boolSC = false; 
opts.boolSA = false; 
opts.boolSN = false;

% Prep the mz vectors
numF = size(data,2);
mzPeaks = cell(numF,1);
for n = 1:numF
    mzPeaks{n,1} = [data(n).mz' nanmean(data(n).sp,1)'];
end

% Determine the cmz vector in the traditional way
cmz = mspmatch(mzPeaks,...
    'estimationMethod',opts.estimationMethod,...
    'mzRes',opts.mzRes,...
    'display',opts.display);

% Now run the samplealign2 function to match A to B...
for n = 1:numF
    
    [j,k] = samplealign2(...
        cmz,...
        mzPeaks{n,1},...
        'Band',opts.handBand,...
        'Gap',opts.handGap,...
        'Distance',opts.handDist,...
        'Quantile',[],...
        'SHOWCONSTRAINTS',opts.boolSC,...
        'SHOWNETWORK',opts.boolSN,...
        'SHOWALIGNMENT',opts.boolSA);

    % Make the matrices to finish...
    newSP = zeros(size(data(n).sp,1),size(cmz,1));
    newSP(:,j) = data(n).sp(:,k);
    
    % Save to the structure
    data(n).al = newSP;

end

cmz = cmz(:,1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

