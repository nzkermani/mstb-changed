function dpnDataExtract(~,~,fig)
%dpnDataExtract - get the annotated pixels from the data

% Guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end
if ~isfield(dpn,'anno')
    disp('No annotations');
    return
end

% Extract the annotated pixels from the data...
switch dpn.mode
    case 'dual'
        [mask2,histID,isInt1,isInt2,pixID] = dpnAnnotationExtract(dpn);
    case 'single'
        [mask2,histID,pixID] = desiAnnotationExtract(dpn);
        isInt1 = [];
        isInt2 = [];
end

% Extract the MS data...
szPos = size(dpn.d1.sp);
spPos = reshape(dpn.d1.sp,[szPos(1)*szPos(2) szPos(3)]);
spPos = spPos(mask2 > 0,:);

% Only do if dual mode
switch dpn.mode
    case 'dual'
        szNeg = size(dpn.d2.sp);
        spNeg = reshape(dpn.d2.sp,[szNeg(1)*szNeg(2) szNeg(3)]);
        spNeg = spNeg(mask2 > 0,:);
    otherwise
        szNeg = [];
        spNeg = [];
end        

% Generate a cell vector of histological ID / class, which is based on the
% text annotations in dpn.anno
numID = mask2(mask2 > 0);
histID = histID(mask2 > 0,:);
numP = numel(numID);

if strcmp(dpn.mode,'dual')
    isInt1 = isInt1(mask2 > 0,:);
    isInt2 = isInt2(mask2 > 0,:);
end

% Need a vector for pixel ID to serve as a proxy for sampleID, as leave one
% out cv doesn't work at the moment

% Now this is where we would consider launching the stats toolbox for some
% data analytical fun
meta1.pixID = pixID(mask2 > 0,:);
meta1.histID = histID;

% All this stuff is only done in dual mode
if strcmp(dpn.mode,'dual')
    meta1.isInt = isInt1;
    meta2.pixID = pixID(mask2 > 0,:);
    meta2.histID = histID;
    meta2.isInt = isInt2;
end

if strcmp(dpn.mode,'dual')
    choice = questdlg('Which ion mode?',...
        'Ion Mode',...
        'Positive','Negative','Negative');
    
    switch choice
        case 'Positive'
            launchStats(dpn.d1.mz,spPos,meta1);
            
        case 'Negative'
            launchStats(dpn.d2.mz,spNeg,meta2);
    end
    
else
    % Single mode analysis
    launchStats(dpn.d1.mz,spPos,meta1);
end

end

