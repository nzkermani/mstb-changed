function desiReadRawCompare
%desiReadRawCompare - dump new bits into existing files for comparative
%purposes

% Folder locations...
oldF = 'E:\Data\Olivia\Neg Reproc Test\Anno\';
newF = 'E:\Data\Olivia\Neg Reproc Test\KV-TB\';
%savF = 'E:\Data\Olivia\Neg Reproc Test\Rep-TB\';
copF = 'Z:\Lab Data\Bioinformatics\Test\WatersRaw\';

% Find file types in old and new folders
fOld = fileFinderAll(oldF,'mat',true);
fNew = fileFinderAll(newF,'mat',true);

% Now perform
for n = 1:size(fNew,1)
    
    % Find matching file name
    fx = ~cellfun(@isempty,strfind(fOld(:,2),fNew{n,2}(1:20)));
    f1 = fOld{fx,2};
    f2 = fNew{n,2};
    
    % Open
    tmp1 = open([oldF f1]);
    tmp2 = open([newF f2]);
    
    dpn = tmp1.dpn
    dpn.d1.mz = tmp2.dpn.d1.mz;
    dpn.d1.sp = tmp2.dpn.d1.sp;
    
    % Save...
    %f3 = [savF f1];
    %save(f3,'dpn');
    
    f4 = [copF f1];
    save(f4,'dpn');
    
    
end



end

