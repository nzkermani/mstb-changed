function stats
%stats - simple function for PCA and multivariate analysis of spectral data
%with limited options for normalisation and stuff like that

% Define a default path
defP = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/';
if ~exist(defP,'dir')
    defP = pwd;
end
if strcmp(defP(end),filesep)
    defP = [defP filesep];
end

% Find existing windows
if strcmp(getUser,'jmckenzi')
    f0 = findobj('Tag','stats','Type','figure');
    delete(f0);
end

% Draw the window to start with
[fig] = statsDraw;

% Set the callbacks
statsCallbacks([],[],fig,defP);


end

