function desiNazaninIMZML(src,event,fig,file,defP,opts)
% desiNazaninIMZML - calling function for the .

% Run the main function, to get the matrix and mz vector
%[X1,MZ1,bin] = desiNazaninWorkflow(file.dir,file.nam,opts);

% Now save the results into the structure as required.
dpn.file = file;
dpn.opts = opts;

dpn.d1.mz = MZ1;
dpn.d1.sp = X1;

dpn.fig = fig;
dpn.defP = defP;

dpn.mode = 'single';

%delete(wb);

% Add the guidata...
guidata(fig.fig,dpn);

% Update the images...
dpnUpdateMS([],[],fig,'force');

end

