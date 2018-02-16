function [devs,files,timeStamp] = imzmlRecal( fold )
%imzmlRecal - import, recalibrate and export new imzML files

% Get imzML files
allF = fileFinderAll(fold,'imzML',true);
numF = size(allF,1);

% Waitbar for feedback
wb = waitbar(0,'imzML recalibration');

opdir = '/Volumes/JSM/DB/Colorectal DESI Recalibrated/';

devs = NaN(6,numF);
files = cell(numF,2);
timeStamp = zeros(numF,1);

% Loop through
for n = 1:numF
    
    % Skip files that already exist...
    tmpSave = [opdir allF{n,2}(1:end-6) '-RECAL' '.ibd'];
    if exist(tmpSave,'file')
        continue;
    else
        %disp('sdfgh');
        %continue;
    end
    
    try
        
        % File name
        imz = [allF{n,1} filesep allF{n,2}];
        disp(imz);
        files{n,1} = allF{n,2};
        files{n,2} = allF{n,1};
        
        % Datestamp
        timeStamp(n,1) = imzmlDate(imz);

        % Import file
        [data] = imzmlRawExtract(imz);

        % Perform recalibration
        %[recal] = imageRecal(data);
        [devs(:,n)] = globalRecal(data,true);
        %[recal] = globalRecal(data);
        
        % In the instance of one specific file, we just trim out the bottom
        % row and be done with it
%         if strcmp(allF{n,2},'LNTO31_16_1.imzML') || strcmp(allF{n,2},'LNTO32_17_2.imzML')
%             recal = recal(1:end-1,:);
%         end

        % Output new imzml file
        %imzmlWrite(imz,recal,opdir);

    catch err
        err
        disp(['FAIL --- ' imz]);
    
    end
        
    waitbar(n/numF,wb,allF{n,2});
    
end

delete(wb);

end

