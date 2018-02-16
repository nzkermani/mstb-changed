function [ info ] = endometrialDates
%endometrialDates - extract acquisition dates from the Z drive

fold = ['/Volumes/Data/Lab Data/Endometrial/Olivia/'...
    'Olivia DESI Raw data/1 scan per sec/'];

% Get the folders in here...
flist = folderList(fold);
numF = size(flist,1);

% Info stored here
info = cell(numF,7);


% Loop through
for n = 1:numF
    
    % Save in info
    info{n,1} = flist{n,1};
    
    % Full folder path
    ff = [fold flist{n,1} filesep];
    
    % Find raw files within
    rawl = fileFinderAll(ff,'raw',true);
    
    % Let's just use the first for simplicity as we assume all acquired on
    % the same day (2 are lockmassed versions anyway)
    tc = datenum(rawl{1,2}(1:8),'yyyymmdd');
    
    info{n,2} = str2double(rawl{1,2}(3:8));
    info{n,4} = str2double(datestr(tc,'mm'));
    info{n,3} = datestr(tc,'mmmm');
    info{n,5} = str2double(datestr(tc,'dd'));
    info{n,6} = str2double(datestr(tc,'yyyy'));
       
    info{n,7} = tc;
    
    info(n,:)
   
    
end
    


end

