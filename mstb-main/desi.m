function desi
%desi - this is the main function for doing DESI data processing

defP = '/Volumes/JSM/DB/';
if ~exist(defP,'dir')
    defP = pwd;
end
if ~strcmp(defP(end),filesep)
    defP = [defP filesep];
end

% Do some simple things related to the imzML converter here...
javaclasspath('packages/imzMLConverter/imzMLConverter.jar');

% Now draw the window
[fig] = desiDraw;

% Set the function callbacks
desiCallbacks([],[],fig,defP);


end

