function [mz,logDS,opts,allChk] = xxxNormTranScal(man,mz,sp,meta)
% xxxNormTranScal - Perform the generic norm/tran/scal procedures
%
% man - the structure to the side menu
% mz  - mz vector of the data
% sp  - the actual data in NON-IMAGE format
%
% xxx becomes the default prefix for all 'common' files that are relevant
% to both desi and desiPosNeg. Going to get confusing.

if nargin == 3
    meta = [];
end

% Need to get the normalisation method from the menu
val2 = get(man.norm,'Value');
lab = get(man.norm,'String');
normMethod = lab{val2};

% What if the normMethod is 'Internal Standard'?
if strcmp(normMethod,'Internal Standard')
    [intStdScFac] = doIntStd(mz,sp);
    
    % If nothing good happened, then just quit
    if isempty(intStdScFac)
        mz = [];
        logDS = [];
        opts = [];
        allChk = [];
        warning('No normalisation performed');
        return
    end    
end

% May need to do filtration based on m/z and RT, rather than just the m/z
% for which this function was originally designed
[mz,allChk,mzLo,mzHi,rtLo,rtHi] = getVariableRange(man,mz);

% Ensure that everything is 0 or above
sp(sp < 0) = 0;

% Perform post-hoc normalisation here...
preNormed = sp(:,allChk);
switch normMethod
    case 'Internal Standard'
        reNormed = bsxfun(@rdivide,preNormed,intStdScFac);
        
    case 'PQN-Year'
        warning('This is not recommended');
        [reNormed] = pqnYearBAD(preNormed,normMethod,meta);
        
    case 'PQN-Global'
        [reNormed,scaleFac] = jsmNormalise(preNormed,'PQN-Mean',0,0,[]);
        
        meta2 = meta;
        meta2.scaleFac = scaleFac;
        assignin('base','meta',meta2);        
        
    otherwise
        [reNormed,~] = jsmNormalise(preNormed,normMethod,0,0,[]);
end

% Apply constant scaling factor
reNormed = 1e4 * reNormed / nanmax(reNormed(:));

% Log offset determination
[logDS,doLog,logOS] = doLogFunction(man,reNormed);

% Now we need to consider the scaling technique to be applied here
[logDS,scale] = doScale(man,logDS);

% Collate the things into a structure and return it
opts.mzRange = [mzLo mzHi];
if isfield(man,'rtL')
    opts.rtRange = [rtLo rtHi];
else
    opts.rtRange = [];
end
opts.normMethod = normMethod;
opts.doLog = doLog == 2;
if opts.doLog
    opts.logOS = logOS;
end
opts.scale = scale;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [reNormed] = pqnYearBAD(preNormed,normMethod,meta)
% Perform PQN normalisation by default, and also subtract the
% year's mean intensity for each variable. This is dubious batch
% correction...
reNormed = preNormed;

% Ensure that we have the date
if ~isfield(meta,'date')
    error('No date can be found');
end

% Convect to year
year = datestr2num(meta.date,'yyyy');
[unq,~,ind] = unique(year);
for n = 1:numel(unq)
    
    % Indices of this year
    fx = ind == n;
    
    % Mean of this year
    my = nanmean(reNormed(fx,:),1);
    
    % Subtract the mean
    reNormed(fx,:) = bsxfun(@minus,reNormed(fx,:),my);
    
end

% Set intensities below zero to be zero
reNormed(reNormed < 0) = 0;

[reNormed,~] = jsmNormalise(preNormed,normMethod,0,0,[]);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [intStdScFac] = doIntStd(mz,sp)

% Default...
intStdScFac = [];

% Draw a menu containing a list of potential peaks...
isMZ = {'169.1020';'173.0540';'184.0661';'503.1627';'885.5498';'OTHER'};

[dlgS,dlgV] = listdlg('PromptString','Select m/z',...
    'SelectionMode','single',...
    'ListString',isMZ,...
    'InitialValue',5);

if dlgV ~= 1    
    return
else
    other = str2double(isMZ{dlgS});
end

% If other was selected, have another box...
if dlgS == numel(isMZ)
    other = inputdlg('m/z','m/z',1,{''});
    if isempty(other)
        return
    else
        % Can we convert the m/z value to a number?
        other = str2double(other{1});
    end
end

% Now we need to find the m/z
if ~isnumeric(other)    
    return
end

% Find the ion within 5 ppm
mzFx = mzFind(mz,other,5);
if sum(mzFx) == 0
    [~,mzFx] = min(abs(mz-other));
end

% This is the scaling factor
intStdScFac = sp(:,mzFx);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,allChk,mzLo,mzHi,rtLo,rtHi] = getVariableRange(man,mz)
% Determine the variable range for analysis

mzLo = str2num(get(man.mzL,'String')); %#ok<*ST2NM>
mzHi = str2num(get(man.mzH,'String'));

if size(mz,2) == 2 && isfield(man,'rtL') && isfield(man,'rtH')
    rtLo = str2num(get(man.rtL,'String'));
    rtHi = str2num(get(man.rtH,'String'));
    mzChk = mz(:,1) >= mzLo & mz(:,1) <= mzHi;
    rtChk = mz(:,2) >= rtLo & mz(:,2) <= rtHi;    
    allChk = mzChk & rtChk;    
    mz = mz(allChk,:);    
else
    rtLo = [];
    rtHi = [];
    allChk = mz >= mzLo & mz <= mzHi;
    mz = mz(allChk);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [logDS,doLog,logOS] = doLogFunction(man,reNormed)
% Perform log transformation

doLog = get(man.log,'Value');
if doLog == 2
    logOS = nanmedian(reNormed(reNormed > 0));
    minOS = min(reNormed(:));
    if minOS+logOS < 1
        logOS = 1;
    end
    logDS = log(reNormed + logOS);
else
    logDS = reNormed;
    logOS = [];
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [logDS,scale] = doScale(man,logDS)
% Perform scaling

sc1 = get(man.scale,'Value');
lab = get(man.scale,'String');
scale = lab{sc1};
switch scale
    case 'UV'
        % Apply UV scaling to give mean 0 and std 1
        logDS = zscore(logDS);
    otherwise
        % Do nothing
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%