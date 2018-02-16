function desiWatersAnalyte(fig,file,defP)
%desiWaters - function to process the Waters Raw files and put the output
%as should be expected for the toolbox

% Do the import / processing function here
[MZ,X] = watersAnalyteFile([file.dir file.nam],[]);

% Update the guidata structure and draw a few things here
dpn.file = file;
dpn.opts = [];

dpn.d1.mz = MZ;
dpn.d1.sp = X;

dpn.fig = fig;
dpn.defP = defP;

dpn.mode = 'single';

% Add the guidata...
guidata(fig.fig,dpn);

% Update the images...
dpnUpdateMS([],[],fig,'force');


end

