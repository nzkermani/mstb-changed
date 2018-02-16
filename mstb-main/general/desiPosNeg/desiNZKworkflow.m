function desiNZKworkflow(~,~,fig,file,~,opts)
%desiNZKworkflow - use functions developed by Nazanin to import a centroid
% mode imzML file.  The options are specified in the desi toolbox. Input
% arguments are:
% src   | callback control; omitted
% event | callback control; omitted
% fig   | structure containing handles to 'desi' gui
% file  | structure with .dir and .nam, containing path to imzML file
% defP  | default path to be saved with file; omitted
% opts  | desi processing options (see desiFileProcess for more info)

% Have a waitbar
wb = waitbar(0.25,'Importing imzML file');

%%%%%%%%%%%%%%%%%%%%% EDIT FUNCTIONS BELOW THIS LINE %%%%%%%%%%%%%%%%%%%%%

% Run the main function, to get the matrix and mz vector
[bin,dims] = nzkPeakMatching3MethodsV2(...
    file.dir,...
    file.nam,...
    opts.ppmRes,...
    opts.mzRange,...
    opts.nzkBL1,...
    opts.nzkSplCor,...
    opts.nzkMethod);

% Update waitbar
waitbar(0.5,wb);

% Filter out those with low frequencies
[mz,sp] = b2iFilter(dims,bin,opts.mzFrac,'withinImage');

% Now remove those which appear mostly on the background rather than the
% actual tissue object. We use the summed ion image over m/z 600-1000
mask = mz > 600 & mz < 1000;
if sum(mask) == 0
    mask = mz > 0;
end
tmpImg = nansum(sp(:,:,mask),3);
[tobg,~] = dpnTOBG(tmpImg,[],[]);
[chk] = nzkRemoveBackgroundIons(sp,tobg,tmpImg);
mz = mz(chk);
sp = sp(:,:,chk);

%%%%%%%%%%%%%%%%%%%%% EDIT FUNCTIONS ABOVE THIS LINE %%%%%%%%%%%%%%%%%%%%%

% Now save the results into the structure as required.
dpn.file = file;
dpn.opts = opts;
dpn.d1.mz = mz';
dpn.d1.sp = sp;
dpn.fig = fig;
dpn.defP = file.dir;
dpn.date = imzmlDate([file.dir file.nam]);
dpn.mode = 'single';

% Delete the waitbar
delete(wb);

% Add the guidata...
guidata(fig.fig,dpn);

% Update the images...
dpnUpdateMS([],[],fig,'force');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
