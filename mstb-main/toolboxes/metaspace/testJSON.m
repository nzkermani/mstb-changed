%% FUNCTION TO TEST imzmlJSON.m

% Specify file path / names
fPath = '/Volumes/JSM/DB/Lymph DESI Recalibrated/';

% Get all of the imzML files from this folder
allF = fileFinderAll(fPath,'imzML',true);
numF = size(allF,1);

% Loop
for n = 1:numF
    
    % Strip out the imzML part...
    fName = allF{n,2}(1:end-6);
    
    imzmlJSON(fPath,fName,...
        'Organism Part','Lymphnode',...
        'Polarity','negative',...
        'Analyser','orbitrap',...
        'Extra','MTBLS385; recalibrated imzML file');
    
end