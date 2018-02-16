function [ fails ] = replaceAnno2Proc(match)
%replaceAnno2Proc - use the list to replace old data with new data. And
%then save to a new location... This works in conjunction with the
%matchAnno2Proc function

% Where are we to save the files?
newSave = 'E:\Data\Olivia\Updated\';

% How many files have been matched?
numF = size(match,1);

flag = false(numF,1);

% Loop through the files
for n = 1:numF
    
    % Where should we save the new file?
    sl = strfind(match{n,1},filesep);
    subFold = match{n,1}(sl(end)+1:end);
    
    % Full new file save name
    newN = [newSave subFold filesep match{n,4}(1:end-4) '-Upd.mat'];
    
    % Check that is doesn't already exist
    if exist(newN,'file')
        disp('Exists already');
        continue;
    end
    
    disp(newN);

    % So what we actually now need to do is load the old annotated matlab
    % file, replace the data and then save in the new location
    dpnAnno = open([match{n,1} filesep match{n,2}]);
    dpnRepr = open([match{n,3} filesep match{n,4}]);
    
    % Check file sizes
    szA = size(dpnAnno.dpn.d1.sp);
    szR = size(dpnRepr.dpn.d1.sp);
    
    if szA(1) == szR(1) && szA(2) == szR(2)
    else
        disp('File size mismatch');
        
        figure;
        subplot(1,2,1);
        imagesc(nansum(dpnAnno.dpn.d1.sp,3));
        title('ORIGINAL');
        subplot(1,2,2);
        imagesc(nansum(dpnRepr.dpn.d1.sp,3));
        title('REPROCESSED');
        
        flag(n,1) = true;
        continue;
    end
        
    
    % This is the file structure
    dpn = dpnAnno.dpn;
    
    dpn.file = dpnRepr.dpn.file;
    dpn.opts = dpnRepr.dpn.opts;
    
    dpn.d1.mz = dpnRepr.dpn.d1.mz;
    dpn.d1.sp = dpnRepr.dpn.d1.sp;
    dpn.d1.numPoints = dpnRepr.dpn.d1.numPoints;
    
    % Save the file...
    if ~exist([newSave subFold],'dir')
        mkdir([newSave subFold]);
    end
    save(newN,'dpn');
    
end

fails = match(flag,:);


end

