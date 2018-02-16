function desiCheck
%% desi toolbox compatibility

% Need a script to determine if your computer has all the requirements to
% run the toolbox.  This is likely to get a little harder when trying to do
% it for Windows and MassLynx compatibility.

%% Perform all of the tests here
path = [pwd filesep 'toolboxes' filesep 'desi' filesep];
name = 'desi.m';

% File exists?
if ~exist([pwd filesep name],'file')
    flFile = false;
else
    flFile = true;
end

% What about mex compilation?
try
    mex([path 'yprime.c']);
    clc;
    flMex = true;
    if ispc
        delete('yprime.mexw64');
    else
        delete('yprime.mexmaci64');
        disp('need to delete file');
    end
    
catch
    flMex = false;
end

% Suitable matlab version
mlv = ver('matlab');
switch mlv.Version
    case {'8.3','9.0'}
        flVer = true;
    otherwise
        flVer = false;
end

% Toolbox list
tbl = ver;
tbList = {tbl.Name}';

% Toolbox lists...
req = {'Bioinformatics';...
    'Image Processing';...
    'Curve Fitting';...
    'Parallel Computing';...
    'Signal Processing';...
    'Statistics and Machine Learning';...
    'Wavelet'};
numT = numel(req);
tbYes = false(numT,1);

% Loop through
for n = 1:numT    
    tbYes(n,1) = any(strcmpi(tbList,[req{n,1} ' Toolbox']));    
end

% % % % Now what about the bioinformatics toolbox?
% % % if any(strcmp(tbList,'Bioinformatics Toolbox'))
% % %     flBio = true;
% % % else
% % %     flBio = false;
% % % end
% % % 
% % % % Image processing toolbox?
% % % if any(strcmp(tbList,'Image Processing Toolbox'))
% % %     flImage = true;
% % % else
% % %     flImage = false;
% % % end
% % % 
% % % % Stats toolbox
% % % 


%% Write a pretty little message to screen to make it look official
clc;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% desiCheck');
disp('%');
disp('% Checking compatibility');
disp('%');

if flFile
    disp(['% Pass' char(9) char(9) 'File check']);
else
    disp(['* Fail' char(9) char(9) 'File check']);
end

if flMex
    disp(['% Pass' char(9) char(9) 'Mex compile']);
else
    disp(['* Fail' char(9) char(9) 'Mex compile']);
end

if flVer
    disp(['% Pass' char(9) char(9) 'Matlab version']);
else
    disp(['* !!!!' char(9) char(9) 'Matlab version warning']);
    disp(['*' char(9) char(9) char(9) char(9) char(9) char(9) 'Only 2014a/2016a have been tested']);
end

% For each of the toolboxes, tell me something
for n = 1:numT
    if tbYes(n,1)
        disp(['% Pass' char(9) char(9) 'Toolbox: ' req{n,1}]);
    else
        disp(['* Fail' char(9) char(9) 'Toolbox: ' req{n,1}]);
    end
end

disp('%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

end
