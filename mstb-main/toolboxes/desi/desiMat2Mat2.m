function [ status ] = desiMat2Mat2(files)
%desiMat2Mat - import reprocessed spectral data from one bunch on .mat
%files into another batch. For good measure, we'll save the results
%elsewhere to prevent files being corrupted.

% Root save location
rootSave = 'C:\Users\jmckenzi\Desktop\OliviaREPROC\';

% How many files?
numF = size(files,1);

% Status
status = cell(numF,3);

% Loop through each file
for n = 1:numF
    
    % Save orig file names...
    status{n,1} = files{n,1};
    status{n,2} = files{n,2};

    % Skip missed out files
    if isempty(files{n,1})
        status{n,4} = 'Absent file';
        continue;
    end
    
    % Load the ORIGINAL file with the annotations...
    ofn = [files{n,1} filesep files{n,2}];
    orig = open(ofn);
    
    % Load the REPROC file
    nfn = [files{n,3} filesep files{n,4}];
    new = open(nfn);
    
    % Simple checks to ensure same size of matrices...
    szO = size(orig.dpn.d1.sp);
    szN = size(new.dpn.d1.sp);
    
    if szO(1) == szN(1) & szO(2) == szN(2)
        
        orig.dpn.d1.mz = new.dpn.d1.mz;
        orig.dpn.d1.sp = new.dpn.d1.sp;
        
        orig.dpn.file = new.dpn.file;
        
        % Save the file...
        newFN = [rootSave files{n,2}];
        dpn = orig.dpn;
        save(newFN,'dpn');
        
        % Update status
        status{n,3} = newFN;
        status{n,4} = 'Pass';
        
    else
        
        % Can't really do anything
        status{n,4} = 'Incorrect sizes';        
    end
        
    
end


end

