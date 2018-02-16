function mzAssignList(mz,spec)
%mzAssignList - given a list of mz values and a reference spectrum of
%intensities, output the assignments...

isTest = false;

% Need to get some of the input parameters - maybe just a menu box to make
% life easier for people who can't type (mztol, threshold, dbPath)
opts.mztol      = {'5'};
opts.threshold  = {'0'};
opts.dbpath     = {'<default>'};
opts.save       = {[pwd filesep 'Annotations-' datestr(now,'yymmdd-HHMMSS')]};

if ~isTest
    pro = {'Annotation tolerance / ppm',...
        'Intensity threshold',...
        'Database',...
        'Save name'};
    def = [opts.mztol,opts.threshold,opts.dbpath,opts.save];
    
    resp = inputdlg(pro,'Provide',repmat([1 100],numel(def),1),def);
    
    opts.mztol = resp{1};
    opts.threshold = resp{2};
    opts.dbpath = resp{3};
    opts.save = resp{4};
    
    % Error checking
    if ~strcmp(opts.dbpath,'<default>') && ~exist(opts.dbpath{1},'file')
        error('Suggested DB cannot be found. Check spelling.');
    end
end

% Apply a threshold to reduce the number of variables needing to be
% annotated
thr = spec >= str2num(opts.threshold); %#ok<*ST2NM>
mz = mz(thr);
%spec = spec(thr);

% Check that there are still some things left...
if isempty(mz)
    error('Threshold is too high!');
end

% Run the function to do the annotation.
switch opts.dbpath
    case '<default>'
        [ass] = peakAnnotationBatch(mz,str2num(opts.mztol));
    otherwise
        [ass] = peakAnnotationBatch(mz,str2num(opts.mztol),opts.dbpath{1});
end

% Dump the data into a text file
[~] = annotation3text(ass,opts.save,mz);

end

