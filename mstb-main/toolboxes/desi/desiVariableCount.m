function [ fInfo,numV ] = desiVariableCount(fold)
%desiVariableCount - determine the number of variables and perhaps other
%information in a bunch of processed mat files for the toolbox

% Which folder?
if isempty(fold)
    fold = '/Volumes/JSM/DB/PBC/PBC Liver/';
end

% Get all of the files within
allF = fileFinderAll(fold,'mat',true);
numF = size(allF,1);

% Save information
fInfo = cell(numF,3);
numV = zeros(numF,1);
tStamp = zeros(numF,1);

% Loop through...
for n = 1:numF
    
    % Open file
    tmp = open([allF{n,1} filesep allF{n,2}]);
    if ~isfield(tmp,'dpn')
        continue;
    end
    
    
    
    % How many variables
    numV(n,1) = size(tmp.dpn.d1.sp,3);
    
    % What is the time stamp?
    
    % File information
    fInfo{n,1} = allF{n,2};
    [~,fInfo{n,2}] = previousFolder(allF{n,1});
    
    disp(int2str(n));
    
end


end

