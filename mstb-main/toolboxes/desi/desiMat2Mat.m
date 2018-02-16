function [ match ] = desiMat2Mat
%desiMat2Mat - import reprocessed spectral data from one bunch on .mat
%files into another batch. For good measure, we'll save the results
%elsewhere to prevent files being corrupted.

% Define folder locations...
orig = 'Z:\Data\Olivia\For_James\';
reproc = 'C:\Users\jmckenzi\Desktop\OliviaALL\';
resave = 'C:\Users\jmckenzi\Desktop\OliviaREPROC\';

% Folder containing other as yet-unprocessed files
other = 'D:\Olivia DESI Batch 2 - 1ss\';
otherFiles = fileFinderAll(other,'raw',true);

% Find all the original files...
origFiles = fileFinderAll(orig,'mat',true);
numO = size(origFiles,1);

% Find the names of the new files - no doubt these won't match perfectly
repFiles = fileFinderAll(reproc,'mat',true);
numR = size(repFiles,1);

% A cell to store matching file names
match = cell(numO,4);

% Let's loop trhough the original file names and find a match
for n = 1:numO
    
    disp(origFiles{n,2});
    
    
    % Separate text according to the three underscored entries
    und = strfind(origFiles{n,2},'_');
    id1 = origFiles{n,2}(1:und(1)-1);
    id2 = origFiles{n,2}(und(1)+1:und(2)-1);
    id3 = origFiles{n,2}(und(2)+1:und(3)-1);
    
    % Can we find this pattern in the repFiles?
    fx1 = ~cellfun(@isempty,strfind(repFiles(:,2),id1));
    fx2 = ~cellfun(@isempty,strfind(repFiles(:,2),id2));
    fx3 = ~cellfun(@isempty,strfind(repFiles(:,2),id3));
    
    % Combine
    fy = fx1 + fx2 + fx3;
    
    % We are looking for a value of 3 to satuisfy the equation
    pick = find(fy == 3);
    if numel(pick) == 1
        match{n,1} = origFiles{n,1};
        match{n,2} = origFiles{n,2};
        match{n,3} = repFiles{pick,1};
        match{n,4} = repFiles{pick,2};
        disp('Perfect match');        
        
    elseif numel(pick) == 0
        
%         % Can we open the original file and read in the raw file name?
%         fullname = [origFiles{n,1} filesep origFiles{n,2}];
%         dpn = open(fullname);
%         
%         % Can we identify this file on Olivia's hard disk?
%         ff = strcmp(otherFiles(:,2),dpn.dpn.file.nam);
%         ff = find(ff);
%         if numel(ff) == 1
%             match{n,1} = origFiles{n,2};
%             match{n,3} = otherFiles{ff,2};            
%         end
        disp('None found');
        
        
        
        
    elseif numel(pick) > 0
        disp('Too many found');
        
        match{n,1} = origFiles{n,1};
        match{n,2} = origFiles{n,2};
        match{n,3} = repFiles(pick,1);
        match{n,4} = repFiles(pick,2);

    end
    
    
    
    
    
end

end

