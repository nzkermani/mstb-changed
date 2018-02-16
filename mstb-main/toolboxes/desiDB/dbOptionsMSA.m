function dbOptionsMSA(~,~,fig,pan)
%dbOptionsMSA - determine all of the options selected from the MSA, and
%then launch into the main function.

% This is the tag for working out which MAT files are which
tag = 'msa';

% Determine the files to be analysed
[files] = getFiles(fig);

% Check that we actually have selected 2 files
if size(files,1) < 2
    eh = errordlg('Need 2 or more files to continue');
    pause(1.5);
    delete(eh);
    return
end

% Check that they actually exist
ext = false(size(files,1),1);
for n = 1:size(files,1)
    if exist([files{n,1} filesep files{n,2}],'file')
        ext(n,1) = true;
    end
end
files = files(ext,:);

% Determine all of the options
 [opts] = getOptions(pan);

% Now run the processing function
[data,files] = dbMSA(files,opts);
data.tag = tag;
data.opts = opts;
data.files = files;

% Save the data as desired, and then we are done and dusted for this
% function.
save([opts.path opts.name '.mat'],'data','tag');

% Also write out a CSV file which can be used for metadata that somebody
% wants to manually add
try
    dbWriteCSV(data.files(:,2),[opts.path opts.name '.csv']);
catch err
    err
    
end

% Now we need to refresh the table
defP = get(fig.fig,'Name');
[allF] = dbFileFind(defP);
[tabDat] = dbFileInfo(allF);
dbTableUpdate(fig,tabDat);
set(fig.fig,'Name',defP);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [files] = getFiles(fig)
% Get the selected files from the table

% Get files from the table
td = get(fig.tab,'Data');

% These are the selected ones
sel = cell2mat(td(:,1)) == 1;
files = td(sel,2:3);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = getOptions(pan)
% Determine all of the relevant options

% m/z resolution
opts.mzRes = str2double(get(pan.mzTol,'String'));

% ppm tolerance
opts.ppmTol = str2double(get(pan.ppmTol,'String'));

% m/z range...
m1 = str2double(get(pan.mzLow,'String'));
m2 = str2double(get(pan.mzHigh,'String'));
opts.mzRange = [m1 m2];

% Polarity
val = get(pan.polarity,'Value');
tmp = get(pan.polarity,'String');
opts.polarity = tmp{val};

% Pixels
val = get(pan.pixel,'Value');
tmp = get(pan.pixel,'String');
opts.pixel = tmp{val};

% Pixel combo
val = get(pan.pixelCombo,'Value');
tmp = get(pan.pixelCombo,'String');
opts.pixelCombo = tmp{val};

% Save path
opts.path = get(pan.path,'String');

% Save name
opts.name = get(pan.name,'String');

% Dates?
opts.dateLow = datenum(get(pan.dateLow,'String'),'ddmmyy');
opts.dateHigh = datenum(get(pan.dateHigh,'String'),'ddmmyy');

% Minimum number of pixels
opts.minPix = str2double(get(pan.minPix,'String'));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%