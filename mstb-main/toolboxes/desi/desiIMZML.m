function desiIMZML(~,~,fig,file,defP,opts)
%desiIMZML - load an imzML file!

% Run the main function, to get the matrix and mz vector
[MZ1,X1,timestamp] = desiFunctionIMZML(file,opts);

% Now save the results into the structure as required.
dpn.file = file;
dpn.opts = opts;

dpn.d1.mz = MZ1;
dpn.d1.sp = X1;

dpn.fig = fig;
dpn.defP = defP;
dpn.date = timestamp;

dpn.mode = 'single';

%delete(wb);

% Add the guidata...
guidata(fig.fig,dpn);

% Update the images...
dpnUpdateMS([],[],fig,'force');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
