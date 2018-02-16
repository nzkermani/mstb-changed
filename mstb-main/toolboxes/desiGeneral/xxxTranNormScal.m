function [mz,logDS] = xxxTranNormScal(man,mz,sp)
% xxxNormTranScal - Perform the generic norm/tran/scal procedures
%
% man - the structure to the side menu
% mz  - mz vector of the data
% sp  - the actual data in NON-IMAGE format
%
% xxx becomes the default prefix for all 'common' files that are relevant
% to both desi and desiPosNeg. Going to get confusing.

% This is the mz range which we want to look at
mzLo = str2num(get(man.mzL,'String')); %#ok<ST2NM>
mzHi = str2num(get(man.mzH,'String')); %#ok<ST2NM>
mzChk = mz >= mzLo & mz <= mzHi;
mz = mz(mzChk);

% Need to get the normalisation method from the menu
val2 = get(man.norm,'Value');
lab = get(man.norm,'String');
normMethod = lab{val2};

% Ensure that everything is 0 or above
sp(sp < 0) = 0;

% Create a prenormed matrix
preNormed = sp(:,mzChk);

% Log offset determination
doLog = get(man.log,'Value');
if doLog == 2
    logOS = nanmedian(preNormed(preNormed > 0));
    logDS = log(preNormed + logOS);
else
    logDS = preNormed;
end

% Remove the offset
logDS = logDS - min(logDS(:));

% Perform post-hoc normalisation here...
[logDS,~] = jsmNormalise(logDS,normMethod,0,0,[]);
logDS     = 1e4 * logDS / nanmax(logDS(:));

% Now we need to consider the scaling technique to be applied here
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
