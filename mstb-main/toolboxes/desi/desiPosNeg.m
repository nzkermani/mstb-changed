function desiPosNeg
%desiPosNeg - this is the main function for doing dual ion mode DESI data
%analysis from a single imzML file.

defP = '/Volumes/JSM/DB/DESI Polarity Switched/IMZML/';
if ~exist(defP,'dir')
    defP = pwd;
end
if ~strcmp(defP(end),filesep)
    defP = [defP filesep];
end

% Do some simple things related to the imzML converter here...
javaclasspath('packages/imzMLConverter/imzMLConverter.jar');

% Now draw the window
[fig] = dpnDraw;

% Set the function callbacks
dpnCallbacks([],[],fig,defP);


end

