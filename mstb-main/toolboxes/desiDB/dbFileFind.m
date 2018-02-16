function [allF] = dbFileFind(path)
%dbFileFind - get all DESI files that are suitable for the new DESI 
% toolbox

% Database location folder - use type to specify which files to be included
type = 'mat';

% What about the processing style? Single mode or pos/neg?
%mode = 'single';

% Let's have a waitbar
wb = waitbar(0,'Finding files');

% Find all .mat files in the folders
allF = fileFinderAll(path,type);
allF(1,:) = [];
numF = size(allF,1);

% Create pass vector to say which files can be included
pvec = true(numF,1);

% % Now determine which of these files is suitable
% for n = 1:numF
%     
%     i = matfile([allF{n,1} filesep allF{n,2}]);
%     fn = fieldnames(i);
%     
%     % Check that file contains the correct structure    
%     if any(strcmp(fn,'dpn'))
%         % This will do for now!
%         pvec(n,1) = true;
%         
%     end
%     
%     waitbar(n/numF,wb);
% end

% Trim the file list
allF = allF(pvec,:);

% Delete the waitbar
delete(wb);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


