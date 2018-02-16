function pdProcess(~,~,fig)
%pdProcess - begin the processing of the individual sections

% Get the guidata
data = guidata(fig.fig);

% Process together or apart?
opt = fig.procOpt.String;
val = fig.procOpt.Value;
opt = opt{val};

% What abot the file name
file = fig.fSelect.UserData;

% What kind of file is it? Depends on the options...
switch lower(file.ext)
    case 'imzml'
        opts = getOptsIMZ;        
    case 'raw'
        opts = getOptsRAW;        
end

% Get the table info
td = fig.segTable.Data;
data.table(:,3:5) = td;

% Tell me how big the data is?
img = fig.ax(2).CData;
sz = size(img);

% Get the section information...
numS = size(data.table,1);
for n = 1:numS
    
    % Section coordinates
    xx = data.table{n,6};
    yy = data.table{n,7};
    
    % Section name
    name = data.table{n,5};
        
    % Do the processing
    switch lower(file.ext)
        case 'imzml';
            
            % Provide the ROI like this - easy
            opts.roi.xx = xx;
            opts.roi.yy = yy;
            
            % Here is the function to process the imzML file
            [mz,sp,timeStamp] = desiFunctionIMZML(file,opts);
      
        case 'raw'
            
            % How does the waters function want the ROI?
            
            % Determine scans according to the region of interest...
            mask = zeros(sz(1),sz(2));
            mask(yy(1):yy(2),xx(1):xx(2)) = 1;
            figure; imagesc(mask);
            mask = reshape(mask,[sz(1)*sz(2) 1]);
            
            % Now we need to get the scanlist, which is basically the xy2D
            % matrix reshapen...
            scl = reshape(data.xy2D,[sz(1)*sz(2) 1]);
            scl = scl(mask == 1);
            
            opts.roi = true;
            opts.roiList = scl;
            
            % Run the extract function
            [mz,sp,~,~,~,opts] = h5waters([file.dir file.nam],opts);
            [timeStamp] = watersTimeStamp([file.dir file.nam]);
            
            
    end
        
    
    % Now we need to format the data into dpn style and then save it using
    % the appropriate file name...
    clear dpn;
    dpn.file = file;
    dpn.opts = opts;
    dpn.d1.mz = mz;
    dpn.d1.sp = sp;
    dpn.fig = [];
    [dpn.defP,~] = previousFolder(file.dir);
    dpn.date = timeStamp;
    dpn.mode = 'single';
    
    % So now we need to save dpn
    save([dpn.defP name '.mat'],'dpn');

    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = getOptsIMZ

opts.method = 'Centroid';

opts.mzFrac   = 0.01;
opts.mzRes    = 0.002;
opts.ppmRes   = 8;
opts.mzRange  = [100 1000];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = getOptsRAW

opts.mzRes      = 0.01;
opts.mzFrac     = 0.01;
opts.peakDetect = 0;
opts.ppmRes     = 10;
opts.mzRange    = [100 1000];
opts.roi        = true;
opts.roiList    = [];
opts.recal      = true;
opts.method     = 'Centroid';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

