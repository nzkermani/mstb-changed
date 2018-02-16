function desiReproc(~,~,fig)
%desiReproc - reprocess the entire file again, using

% guidata and non-standard return
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Determine the raw file that needs to be reprocessed.
% extn = dpn.file.ext;
% switch extn    
%     case 'imzML'        
%         fn = [dpn.file.dir dpn.file.nam];
%     case 'raw'
%         fn = 'fff';
%         
%     case 'h5'
%         fn = [dpn.file.dir];
%         dpn.file.ext = fileExtension(fn);
%         if strcmp(dpn.file.ext(end),filesep)
%             dpn.file.ext = dpn.file.ext(1:end-1);
%         end
%         
%     case 'h5-raw'
%         fn = [dpn.file.dir dpn.file.nam];
%         dpn.file.ext = 'raw';
%     otherwise
%         disp('There is no otherwise');
%         return
% end
%     
% disp('FILE TO BE RE-PROCESSED');
% disp(fn);

% Check that the file exists... If not, then ask for another one.
%if ~exist(fn,'file')
    [f1,f2,flag] = uigetfile({'*.imzML; *.dat; *.DAT; *.mat'},...
        'File Select',dpn.file.dir);
    fn = [f2 f1];
%end

% Ensure that a file was selected and exists...
if isnumeric(f2)
    return;
end

% Get the extension of the selected file
extn = fileExtension(fn);


% Get the processing options
switch lower(extn)
    case {'raw','imzml','dat'}
        % Now we need to determine the options
        dpn.file.dir = f2;
        dpn.file.nam = f1;
        dpn.file.ext = extn;
        
        % Convert if necessary
        if strcmpi(extn,'dat')
            dpn.file = dat2raw(dpn.file);
        end

        % Get the options
        [opts,flag] = desiGetProcOptions(dpn.file);
        
        if ~flag
            return
        end
        
    case 'mat'
        % Then we don't need to do anything here
    otherwise
        % Is there an otherwise?
        disp('Otherwise?');
end

% Reformat the .dat files to the .raw file name instead
switch lower(extn)
    case 'dat'
        fn = f2(1:end-1);
        extn = 'raw';
    otherwise
        disp('no otherwise');        
end

% Now we can reprocess the file using these parameters
switch lower(extn)
       
    case 'raw'
        % This is the basic Waters/RAW processing function
        [MZ,X,~,~,numPoints,opts] = h5waters(fn,opts.raw);
        
        dpn.d1.mz = MZ;
        dpn.d1.sp = X;
        dpn.d1.numPoints = numPoints;
        dpn.opts = opts;
        dpn.file.name = fn;
        
        
    case 'imzml'
        [MZ1,X1] = desiFunctionIMZML(dpn.file,opts.imzml);
        dpn.d1.mz = MZ1;
        dpn.d1.sp = X1;
        dpn.opts = opts.imzml;
        
    case 'mat'
        % Just read in the files from the new one and replace the sp data
        % matrix...
        tmp = open(fn);
        
        % Make sure that the images are the same size!
        szOrig = size(dpn.d1.sp);
        szNew  = size(tmp.dpn.d1.sp);        
        chk = szOrig(1:2) == szNew(1:2);
        if sum(chk) ~= 2
            error('Images not same size');
        end
        
        % Replace the opts, d1.sp, d1.mz, d1.numPoints
        dpn.d1.mz = tmp.dpn.d1.mz;
        dpn.d1.sp = flipud(tmp.dpn.d1.sp);
        
        % numPoints is relatively new
        try
            dpn.d1.numPoints = flipud(tmp.dpn.d1.numPoints);
        catch            
            dpn.d1.numPoints = [];
        end
        
        % And the processing options
        dpn.opts = tmp.dpn.opts;
        
        % Ditch the extra stuff
        clear tmp;
        
    otherwise
        disp('Currently unsupported');
        return
end

% Now that we have the reprocessed data, we should save it to the structure
% and then we have the new data...
guidata(dpn.fig.fig,dpn);

% Update the images...
dpn.d1.img = log2(nansum(dpn.d1.sp,3)+1);
dpnIonImage([],[],dpn.fig.ax.ms1,dpn.d1.img);


end

