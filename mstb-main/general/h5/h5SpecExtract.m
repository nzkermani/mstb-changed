function [data] = h5SpecExtract(path)
% h5SpecExtract - extract data from the h5 files predominantly with the
% reason to look at the raw spectra before they get processed. Use it as an
% option to review the spectra by plotting them.

mzRange = [600 1000];

% Ask the user for a folder in which to find the h5 files...
if isempty(path)
    path = uigetdir('/Volumes/Data/Data/Breast/');
    if isempty(path)
        return
    end
end
disp(path);

% Find all files
[fList] = findFileType(path,'h5');
[fList] = checkFiles(fList);

% Delete some troublesome files
%fList([1 3],:) = [];
%fList = fList(1:6,:);

numF = size(fList,1);

% LEt's have a waitbar
wb = waitbar(0,'Initialising','Name','Meta Process');

% Initialise storage
[data] = initialiseStorage;

% Loop through each file
for n = 1:numF
    
    % Update the waitbar...
    waitbar(n/numF,wb,fList{n,2});
    
    % Prepare file name
    data(n).fP = [fList{n,1} filesep];
    data(n).fN = fList{n,2};
    fName = [data(n).fP data(n).fN];
    disp(data(n).fN);

    % Read in the spectra
    [data(n).mz,data(n).full,data(n).size] = importH5parts(fName,mzRange);
    
    % Get the annotated regions (and MMC predictions)
    [data(n).histIDNames,data(n).annoROI,data(n).probMMC] = getAnnotations(fName);

    % Extract annotated pixels using the ROIs
    [data(n).X,data(n).histIDs,...
        data(n).sampleIDs,data(n).patientIDs] = annotatedExtract(data(n),'annotated');
    
    % Make the full data sparse to save space (or remove it)
    data(n).full = [];%sparse(data(n).full);
    
    % Let's run some QC checks here to ensure that we don't keep empty
    % spectra or unnecessary pixels
    [data(n)] = spectralQC(data(n),20,0.05);
    
    % Maybe a group plot?
    %groupPlot(data(n).mz,data(n).X,data(n).histIDs);
        
end

delete(wb);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fList] = checkFiles(fList)

numF = size(fList,1);
pass = ones(numF,1);

for n = 1:numF
    
    fN = [fList{n,1} filesep fList{n,2}];
    
    try
        h5info(fN,'/mz');
        h5info(fN,'/X');
    catch
        disp(['XXX' char(9) fList{n,2}]);
        pass(n,1) = 0;
    end

end

fList = fList(pass == 1,:);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = initialiseStorage
% A structure to store all the data

data = struct('fP',[],'fN',[],...
    'size',[],...
    'mz',[],...
    'X',[],...
    'annoROI',[],...
    'probMMC',[],...
    'full',[],...
    'sampleID',[],...
    'patientID',[],...
    'inds',[],...
    'sampleIDs',[],...
    'patientIDs',[],...
    'histIDs',[],...
    'histIDNames',[],...
    'histIDNum',[]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,tmp,sz] = importH5parts(fName,mzRange)

% Read in the mz vector
mz = h5read(fName,'/mz');

mask = mz >= min(mzRange) & mz <= max(mzRange);
mz = mz(mask);

% Read in the full data matrix
tmp = h5read(fName,'/X');
tmp = tmp(:,:,mask);

% And unlog it
tmp = exp(tmp) - min(exp(tmp(:)));

% If there is still a non-zero minimum value we need to act
minVal = min(tmp(:));
if minVal > 0
    disp(['NonZERO: ' fName]);
end

% Reshape it to a 2D matrix
sz = size(tmp);
tmp = reshape(tmp,[prod(sz(1:2)) sz(3)]);

% This just rescales the data to put it all on a simliar scale for when it
% gets to the combination stage of the procedure
maxI = max(tmp(:));
tmp = 5000 * tmp / maxI;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [att,roi,pix] = getAnnotations(fName)

% Annotation names
info = h5info(fName,'/tissue_id');    
numAtt = size(info.Attributes,1) - 1; % one is 'loc'
att = cell(numAtt,2);
for r = 1:numAtt    
    att{r,1} = info.Attributes(r).Name;
    att{r,2} = info.Attributes(r).Value;
end

% These are the rectangluar regions
roi = h5read(fName,'/groupPixels');

% These are the classified pixels
pix = h5read(fName,'/pixelProbs');

% Here do the reshaping to put everything in 2 dimensions
sz = size(roi);
roi = reshape(roi,[prod(sz(1:2)) size(roi,3)]);
pix = reshape(pix,[prod(sz(1:2)) size(pix,3)]);
att = att(:,2);

% Remove the background pixels
bgidx = strcmpi(att,'bg') | strcmpi(att,'background');
roi(:,bgidx) = [];
pix(:,bgidx) = [];
att(bgidx,:) = [];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spec,hists,samples,patient] = annotatedExtract(data,type)
% Extract the annotated (predicted) pixels from the full size data

switch lower(type)
    case 'annotated'
        idx = data.annoROI == 1;
    
    case {'predicted','mmc'}
        idx = data.probMMC > 0.9;
end

% idx tells us which pixels were annotated. check that there are no double
% annotated/predicted pixels
chk = sum(idx,2) > 1;
idx(chk,:) = 0;

extr = struct('hist',[],'spec',[],'patient',[],'sample',[]);
numG = size(idx,2);

for n = 1:numG
    
    % These pixels from class n
    fx = idx(:,n) == 1;
    
    extr(n).spec = data.full(fx,:);
    extr(n).hist = repmat(data.histIDNames(n,1),[sum(fx) 1]);
    
    % Additional information for every pixel
    extr(n).sample = repmat({data.fN},[sum(fx) 1]);
    
    extr(n).patient = extr(n).sample;%repmat({data.fN},[sum(fx) 1]);
    
end

% Concatenate it all
hists = vertcat(extr.hist);
spec  = vertcat(extr.spec);
samples = vertcat(extr.sample);
patient = vertcat(extr.patient);

return

[grps,~,inds] = unique(hists);

pcs = 1e4 * bsxfun(@rdivide,spec,nansum(spec,2));
%pcs = jsmNormalise(spec,'pqn-median',0,0);
los = median(pcs(pcs > 0));
pcs = log(pcs + los);

[ll,ss,ee] = princomp(pcs,'econ');
figure; scatter(ss(:,1),ss(:,2),80,inds,'o','filled');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = spectralQC(data,minPeaks,lowPop)
% Run some simple checks to determine if some observations are rubbish or
% if some m/z variables are likely to be noise based on poor consistency
% through the sample

% How many non-zero peaks for each observation?
nzp = sum(data.X > 0,2);
nzpp = nzp / size(data.X,2);

% Remove observations which have 20 or less peaks
ditch = nzp <= minPeaks;
disp(['Deleting ' int2str(sum(ditch)) ' observations']);

% Now delete these observations
data.X(ditch,:) = [];
data.sampleIDs(ditch,:) = [];
data.histIDs(ditch,:) = [];
data.patientIDs(ditch,:) = [];

% Let's run the identifyTails function to determine if variables are tails
% of larger peaks
[tails] = identifyTails(data.mz,data.X,data.histIDs,false);
data.X = data.X(:,~tails);
data.mz = data.mz(~tails);
disp(['Deleting ' int2str(sum(tails)) ' variables']);
disp('------------------------------------------------------------');


% THis function is perhaps best left until after the peaks have been
% matched
%[didx] = msFindLowPop(data.X,data.histIDs,lowPop);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = archive
% Loop through each...
for n = 1:numF
    
    % Update the waitbar...
    waitbar(n/(1+numF),wb,['Preparing data. Sample '...
        int2str(n) '/' int2str(numF)]);
    
    
    if isempty(fName)
        continue;
    end
    
    % Need to move the reading in of data to above this... so that we can
    % do pixel-extraction in the new manner, i.e. to do a bit of a QC check
    % on pixels before we decide to keep them all
    
    % Read in the mz vector
    mz = h5read(fName,'/mz');
    
    % Read in the full data matrix and unlog it
    tmp = h5read(fName,'/X');
    tmp = exp(tmp) - min(exp(tmp(:)));
    
    % If there is still a non-zero minimum value we need to act
    minVal = min(tmp(:));
    if minVal > 0
        disp(['NonZERO: ' fName]);
    end
    
    % Reshape it to a 2D matrix
    sz = size(tmp);
    dbdata(n).size = sz;
    tmp = reshape(tmp,[prod(sz(1:2)) sz(3)]);

    % Is there a limited m/z range?
    opt = getOption(opts,'m/z range');
    if opt.do == 1
        range  = [str2double(opt.defparams{1,1}{1,2})...
            str2double(opt.defparams{1,1}{1,4})];
    else
        range = [1 Inf];
    end    
    mzinds = mz >= range(1) & mz <= range(2);

    % Are we progressing with the MMC predicted, or the annotated only
    % regions?
    opt = getOption(opts,'Pixel source');
    optName = opt.methodnames{opt.selected == 1};
        
    % Noticed a problem with the selection of the annotation boxes
    % if numel(optName) == 2
    %     delete(wb);
    %     error('Problem with the Source of Raw Data options');
    % end
    
    switch optName
        case 'MMC predicted pixels'
            dbdata(n).inds = files(n).mmcVar;
            dbdata(n).histIDNames = files(n).mmcLab;
            
        case 'Annotated pixels'
            dbdata(n).inds = files(n).roiVar;
            dbdata(n).histIDNames = files(n).roiLab;
            
        case 'QC of annotated pixels'
            
            % Here we need to do the QC check on the pixels            
            dbdata(n).inds = files(n).roiVar;
            dbdata(n).histIDNames = files(n).roiLab;
            
            % We need to do it for individual tissue types
            [unqT,~,indT] = unique(dbdata(n).histIDNames);
            numTT = numel(unqT);
            finalKeep = zeros(size(indT));
            for cnt = 1:numTT
                fx = indT == cnt;
                
                % These are the actual indices
                useIdx = dbdata(n).inds(fx);
                
                if strcmpi(unqT{cnt},'background') || strcmpi(unqT{cnt},'bg')
                    [include] = specCheck(tmp(useIdx,mzinds),true);
                else
                    [include] = specCheck(tmp(useIdx,mzinds));
                end
                
                % Add in to the master matrix of kept pixels
                finalKeep(fx) = include;
            end
            
            % Final check to ensure that 885.55 peak is non-zero...
            % if range(1) < 885 & range(2) > 887
            %     [~,pi38] = min(abs(mz - 885.55));
            %     pisum = tmp(dbdata(n).inds,pi38) > 0;
            % 
            %     % Combine the two
            %     finalKeep = finalKeep & pisum;
            % end
            
            % Now trim out the undesired pixels
            finalKeep = finalKeep == 1;
            dbdata(n).inds = dbdata(n).inds(finalKeep,:);
            dbdata(n).histIDNames = dbdata(n).histIDNames(finalKeep,:);
            
    end
    
    
       
    % What about the random selection of pixels according to the original
    % function? This trimmed down the number of pixels by setting inds to
    % zero for undesired pixels
    
    % Get the m/z vector
    dbdata(n).mz    = mz(mzinds);
    
    % Prepare the actual spectral data & determine the sum of pixels
    % LEGACY: dbdata(n).full  = sparse(tmp(inds,mzinds)); 
    tmpDM = sparse(tmp(dbdata(n).inds,mzinds));
    zchk = nansum(full(tmpDM),2) == 0;
    disp(['file ' int2str(n) ', num zero = ' int2str(sum(zchk))]);
    if sum(zchk) > 0
        
        % Trim the data matrix
        tmpDM = tmpDM(~zchk,:);
        dbdata(n).X = tmpDM;
        
        % Trim observation names and indices...
        dbdata(n).inds = dbdata(n).inds(~zchk);
        dbdata(n).histIDNames = dbdata(n).histIDNames(~zchk);
        
    else
        % Do this as normal for files in which no pixels were zero
        dbdata(n).X = tmpDM;        
    end
                
    % Fill in other perhaps useful information here...
    dbdata(n).sampleID  = files(n).fileName;
    dbdata(n).fP        = files(n).filePath;
    dbdata(n).fN        = files(n).fileName;
    
    % Get patient ID - this is Ottmar's addition, which needs to be
    % maintained
    try
        metaInfo = h5info(fName,'/meta');        
        pIDNum = find(ismember({metaInfo.Attributes.Name},'patient_ID'));
        if min(size(pIDNum)) == 0
            dbdata(n).patientID = dbdata(n).sampleID;
        else
            dbdata(n).patientID = metaInfo.Attributes(pIDNum).Value;
        end
    catch
        dbdata(n).patientID = {dbdata(n).sampleID};
        disp('patient_ID was not specified in metadata');
    end
        
    % Stuff about patient ID and sample ID need to be revised with Nima, as
    % LPO-CV relies on knowing the patients from each other...
    if isnumeric(dbdata(n).patientID)
        %dbdata(n).patientIDs = repmat(dbdata(n).patientID,numel(inds),1);
        dbdata(n).patientIDs = cellstr(num2str(repmat(dbdata(n).patientID,numel(dbdata(n).inds),1)));
    else
        dbdata(n).patientIDs = cellstr(repmat(dbdata(n).patientID,numel(dbdata(n).inds),1));
    end
    
    if isnumeric(dbdata(n).sampleID)
        %dbdata(n).sampleIDs  = repmat(dbdata(n).sampleID,numel(inds),1);
        dbdata(n).sampleIDs  = num2cell(repmat(dbdata(n).sampleID,numel(dbdata(n).inds),1));
    else
        dbdata(n).sampleIDs  = cellstr(repmat(dbdata(n).sampleID,numel(dbdata(n).inds),1));
    end
            
end

% Now that we have read everying in, we need to convert the histIDNames to
% a set of consistent numbers across the full range, i.e. rather than a
% list of names, convert these to numbers and a reduced legend...
allLabs = vertcat(dbdata.histIDNames);
unqLabs = unique(allLabs);

for n = 1:numF
    
    newID = zeros(size(dbdata(n).histIDNames,1),1);
    for r = 1:numel(unqLabs)
        fx = strcmp(dbdata(n).histIDNames,unqLabs{r});
        newID(fx,1) = r;        
    end
    dbdata(n).histIDs = dbdata(n).histIDNames;
    dbdata(n).histIDNames = unqLabs;
    dbdata(n).histIDNum = 1:numel(unqLabs);
        
end

delete(wb);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


