function [fails] = mirrorAnno2Proc
%replaceAnno2Proc - use the list to replace old data with new data. And
%then save to a new location... This works in conjunction with the
%matchAnno2Proc function

% Where is the mirror folder located?
fold = 'E:\Data\Olivia\Z Mirror Pos\';
%fold = 'E:\Data\Olivia\Z Mirror\';

% Where to save to on the Z drive
zSave = 'Z:\Lab Data\Endometrial\Olivia\Olivia DESI Raw data\1 scan per sec\';

% We need to get the dir list for this folder, then loop through it...
list = dir(fold);
numF = size(list,1);

% Fail list
fails = cell(numF,2);

% Loop through...
for n = 1:numF
    
    % Is it a proper folder?
    nam = list(n).name;
    if length(nam) < 3 || list(n).isdir == 0
        disp(nam)
        continue;
    end
    
    % So now we can determine the mat files in this folder
    thisF = [fold nam filesep];
    matf = dir(thisF);
    
    % Concat file names
    allF = {matf.name}';
    origFiles = allF;
    
    % Find the interesting files...
    repr = ~cellfun(@isempty,strfind(allF,'Reproc.mat'));
    if sum(repr) == 1
        repF = allF{repr};
        allF(repr) = {'<Taken>'};
    else
        repF = '';
        fails{n,1} = nam;
        fails{n,2} = 'No reproc file';
        continue;
    end
    
    % Can we find the originally annotated file
    orig = ~cellfun(@isempty,strfind(allF,'.mat'));
    orig2 = ~cellfun(@isempty,strfind(allF,'POS')) | ...
        ~cellfun(@isempty,strfind(allF,'pos'));
    orig = orig & orig2;
    if sum(orig) == 1
        origF = allF{orig};
        allF(orig) = {'<Taken>'};
    else
        fails{n,1} = nam;
        fails{n,2} = 'No original annotations';
        origF = '';
    end
    
    disp(['ORIG = ' origF]);
    disp(['NEW  = ' repF]);
    
    % There are two paths available to us. Either the one where we have
    % both an original annotated mat file along with a newly reprocessed
    % file. Or just the newly reprocessed file.
    if length(origF) == 0
        flag = false;
    else
        flag = true;
    end
    
    % First we need to generate a new file name for the brand new mat file,
    % whether it has annotations or not
    %jsmFile = [thisF nam '_Neg_JSM.mat'];
    %jsmFile = ['E:\Data\Olivia\Z Mirror Test\' nam '_Neg_JSM.mat'];
    zSaveFile = [zSave nam filesep nam '_Pos_JSM.mat'];
    
    if ~any(strfind(nam,'T13_'))
        continue;
    end
    disp(zSaveFile);
    
    % Open the newly processed file
    dpnN = open([thisF repF]);
    
    % Can we open the original annotated file (if it exists)
    if flag
        dpn = open([thisF origF]);
        
        % Check file sizes
        szA = size(dpn.dpn.d1.sp);
        szR = size(dpnN.dpn.d1.sp);
        if szA(1) ~= szR(1) || szA(2) ~= szR(2)
            disp('File size mismatch');

            figure;
            subplot(1,2,1);
            imagesc(nansum(dpn.dpn.d1.sp,3));
            title('ORIGINAL');
            subplot(1,2,2);
            imagesc(nansum(dpnN.dpn.d1.sp,3));
            title('REPROCESSED');
            
            fid = fopen([thisF nam '_FILE_SIZE_MISMATCH.txt'],'w');
            fclose(fid);
            
            fails{n,1} = nam;
            fails{n,2} = 'File size mismatch';
            continue;
        end

        % Now replace bits of dpn with bits of dpnN
        dpn.file = dpnN.dpn.file;
        dpn.opts = dpnN.dpn.opts;
    
        dpn.d1.mz = dpnN.dpn.d1.mz;
        dpn.d1.sp = dpnN.dpn.d1.sp;
        dpn.d1.numPoints = dpnN.dpn.d1.numPoints;

        
    else
        % If no annotations, then save dpnN as dpn
        dpn = dpnN;
        clear dpnN;
    end
    
    % Write the new dpn to file
    dpn = dpn.dpn; %#ok<NASGU>
    %save(jsmFile,'dpn');
    
    % Save remotely to the Z drive
    save(zSaveFile,'dpn');
    
    
    
end



end

