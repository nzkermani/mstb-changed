function desiWaters( file,opts,fig,defP )
%desiWaters - function to process the Waters Raw files and put the output
%as should be expected for the toolbox

% Do the import / processing function here
[MZ,X,~,~,numPoints,opts] = h5waters([file.dir file.nam],opts);

% Update the guidata structure and draw a few things here
dpn.file = file;
dpn.opts = opts;

dpn.d1.mz = MZ;
dpn.d1.sp = X;
dpn.d1.numPoints = numPoints;

dpn.fig = fig;
dpn.defP = defP;

dpn.mode = 'single';

% Add the guidata...
guidata(fig.fig,dpn);

% Update the images...
dpnUpdateMS([],[],fig,'force');

% Now show the QC image
desiQCimages(dpn.opts.qc,true);

end

