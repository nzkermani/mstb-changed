function desiH5Load(fig,file,defP)
%desiH5Load - load an existing h5 file into the desi toolbox

% This function actually does the bulk of the work, and can be called
% independently of the GUI
[mz,x,img,opt,allAnno,imzML] = desiH5LoadFunction(file);

% Create a guidata structure, and then save the results, and then display
% them in the right boxes...
dpn.file = file;
dpn.file.imzML = imzML;

% This bit is a little unknown at the moment. Presumably they need to be
% specific to this toolbox...
dpn.opts = [];

% Copy the mz and sp vector/matrices
dpn.d1.mz = mz;
dpn.d1.sp = x;
dpn.d1.img = img;

% Coregistered images...
dpn.opt.coreg = opt.coreg;

% These were provided
dpn.fig = fig;
dpn.defP = defP;

% And the annotations
dpn.anno = allAnno;

% Mandatory for single non PS data
dpn.mode = 'single';

% From here we could do the dpnMatlabLoad function, just skipping the
% reading in functionality part...
dpnMatlabLoad(dpn,dpn.fig);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%