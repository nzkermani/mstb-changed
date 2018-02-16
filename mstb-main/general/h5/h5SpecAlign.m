function [data,x,avg] = h5SpecAlign(data,x)
%h5SpecAlign - perform alignment of the data, form the matrices, and then
%move on to normalisation


% PEAK MATCHING
if ~isfield(data(1),'almz')
    
    maxPS = 8;
    
    % This is where we do some peak matching
    [data,~] = msipmatch(data,...
        'mzRes',0.001,...           % Da
        'maxPeakShift',maxPS,...        % ppm, uses the handle
        'gapPenalty',maxPS,...          % ppm, uses the handle
        'correctionMethod','shortest-path',....
        'estimationMethod','histogram');

    % Diagnostics plot for alignment
    alignmentDiagnostics(data,0.001/10,false,[]);
    
    x = [];
    avg = [];
    return

end

% SINGLE MATRIX CREATION
if isempty(x)
    [x] = createMatrices(data,data(1).almz,'otherwise');
end

% RARE VARIABLE FILTRATION - perhaps before combination of spectral
% average?
[didx] = msFindLowPop(x.X,x.histID,0.025);
x.X = x.X(:,~didx);
x.mz = x.mz(~didx);


% COMBINE SPECTRA
avg.mz = x.mz;
[avg.X,avg.sampleID,...
    avg.histID,idx] = detAvgSpectra(x.X,x.sampleID,x.histID,'mean');%,200);
avg.patientID = x.patientID(idx,:);






end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x] = createMatrices(data,cmz,method,n2a)
% Go from the separate but peak-matching data matrices to a single dataset
% that is ready for analysis

% Create data matrix from structure fields
x.mz        = [];
x.X         = sparse(vertcat(data.align));
x.histID    = vertcat(data.histIDs);
x.patientID = vertcat(data.patientIDs);
x.sampleID  = vertcat(data.sampleIDs);


% Remove empty variables
zeroIndices = sum(x.X,1) == 0;
x.X(:,zeroIndices) = []; 
cmz = cmz(:,1)';
cmz(zeroIndices) = [];
x.mz = cmz;

return

% Now let's integrate pixels by averaging etc...
%opt = getOption(opts,'Pixel selection');
%method = opt.methodnames{opt.selected == 1};
switch lower(method)
    
    case 'n-spectral mean (per sample)'
        %n2a = str2double(opt.defparams{opt.selected == 1}(2));
        [XPeaks,sampleID,patientID,histID] = integTechReplMZ2(XPeaks,...
            sampleID,patientID,histID,n2a,'mean');
        
    case 'n-spectral median (per sample)'
        %n2a = str2double(opt.defparams{opt.selected == 1}(2));
        [XPeaks,sampleID,patientID,histID] = integTechReplMZ2(XPeaks,...
            sampleID,patientID,histID,n2a,'median');
        
    case 'n pixels (per class per file)'
        %n2a = str2double(opt.defparams{opt.selected == 1}(2));
        [ssi] = nPixPerClassPerFile(n2a,sampleID,histID);
        
        % Now we need to decide about the ssi and trim the data
        x.cmz         = cmz;
        x.XPeaks      = XPeaks(ssi,:);
        x.sampleID    = sampleID(ssi,:);
        x.patientID   = patientID(ssi,:);
        x.histID      = histID(ssi,:);
        return
        
    case 'n pixels (per class)'
        %n2a = str2double(opt.defparams{opt.selected == 1}(2));
        [ssi] = nPixPerClass(n2a,sampleID,histID);
                
        % Now we need to decide about the ssi and trim the data
        x.cmz         = cmz;
        x.XPeaks      = XPeaks(ssi,:);
        x.sampleID    = sampleID(ssi,:);
        x.patientID   = patientID(ssi,:);
        x.histID      = histID(ssi,:);
        return

    otherwise
        % Use all of the pixels! So do nothing        
end       
    
% Place into a single structure for the output
x.cmz       = cmz;
x.XPeaks    = XPeaks;
x.sampleID  = sampleID;
x.patientID = patientID;
x.histID    = histID;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ssi] = nPixPerClassPerFile(n2a,sampleID,histID)
% Randomly take n pixels for each class from each file

% For each slide, we determine the unique classes, then select a random
% number n2a of pixels
[unqS,~,indS] = unique(sampleID);
ssi = zeros(size(sampleID,1),numel(unqS));

for s = 1:numel(unqS)

    fx = indS == s;

    % These are the unique histIDs in this sample
    [unqC,~,indC] = unique(histID(fx,:));

    % Now for each of the groups we need to find the number of
    % pixels in each class
    for n = 1:numel(unqC)

        % These are the pixels of class n in sample s
        fy = indC == n;
        if sum(fy) <= n2a
            % Then we do nothing - all pixels are to be selected
            fz = fx & strcmp(histID,unqC{n});
            ssi(:,s) = ssi(:,s) + fz;

        else
            % Then we need to do a subset
            fz = fx & strcmp(histID,unqC{n});
            fz = find(fz);

            % Select a subset
            fa = randperm(numel(fz),n2a);

            % Add into the overall matrix
            ssi(fz(fa),s) = ssi(fz(fa),s) + 1;

        end

    end

end

ssi = sum(ssi,2) > 0;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ssi] = nPixPerClass(n2a,sampleID,histID)
% Randomly take n pixels for each class over all files


% For all slides we determine the unique classes, then take an even number
% of pixels of that class from each of the files in which that class is
% found
[unqC,~,indC] = unique(histID);
ssi = zeros(size(sampleID,1),numel(unqC));

for n = 1:numel(unqC)
    
    fx = indC == n;
    
    % Which files is this found in?
    [unqS,~,indS] = unique(sampleID(fx,:));

    % How many files?
    numF = numel(unqS);
    
    % Find out how many pixels in each file
    ppf = zeros(numF,1);
    for r = 1:numF        
        % Indices of matching histo / sample
        fy = fx & strcmp(sampleID,unqS{r});        
        ppf(r,1) = sum(fy);
    end
    
    % Now decide how many from each of the files
    maxP = round(n2a / numF);
    ppf = min(ppf,repmat(maxP,numF,1));
    
    % Now go and get the actual pixels...
    for r = 1:numF
        fy = fx & strcmp(sampleID,unqS{r});
        fy = find(fy);
        
        if numel(fy) <= ppf(r)
            % Then use all pixels
            ssi(fy,n) = ssi(fy,n) + 1;
        else
            % Then subset
            fa = randperm(numel(fy),ppf(r));
            ssi(fy(fa),n) = ssi(fy(fa),n) + 1;
        end
        
    end 
    
end
    
ssi = sum(ssi,2) > 0;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
