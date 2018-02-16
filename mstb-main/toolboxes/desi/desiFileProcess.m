function desiFileProcess(~,~,fig,man,file,defP)
%desiFileProcess

% Get the options from the processing menu - should be the same for both
% raw and imzML files (up to a point)
opts = getOptions(man);

% What about the various files?
switch opts.file
    
    case 'imzML'
        
        % What about the imzML method
        switch opts.method
            case 'Normal'
                desiIMZML([],[],fig,file,defP,opts);
                
            case 'Nazanin'                
                desiNZKworkflow([],[],fig,file,defP,opts);                
                
        end
        
    case 'RAW'
        
        % What about profile / centroid mode?
        switch opts.method
            
            case 'Profile'
                opts.peakDetect = 1;
                opts.recal      = 0;
            case 'Centroid'
                opts.peakDetect = 0;
                opts.recal      = 1;                
        end
        
        desiWaters(file,opts,fig,defP);
        
    otherwise
        error('There is no otherwise')
end
        
    



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = getOptions(man)

% What is the file type?
tmp = get(man.file,'Value');
val = get(man.file,'String');
opts.file = val{tmp};

% What is the value for method?
tmp = get(man.method,'Value');
val = get(man.method,'String');
opts.method = val{tmp};

opts.mzFrac   = str2double(get(man.mzFrac,'String'));
opts.mzRes    = str2double(get(man.mzRes,'String'));
opts.ppmRes   = str2double(get(man.ppm,'String'));
opts.mzRange  = [str2double(get(man.mzL,'String')) ...
    str2double(get(man.mzH,'String'))];

% Get the ROI value
tmp = get(man.roi,'Value');
val = get(man.roi,'String');
if strcmp(val{tmp},'On')
    opts.roi = true;
    try
        opts.roiList = evalin('base','roi');
    catch
        opts.roiList = [];
    end
else
    opts.roi = false;
    opts.roiList = [];
end

% Get the other options...
if strcmp(opts.method,'Nazanin')
    
    opts.nzkMethod = man.nzkMethod.String{man.nzkMethod.Value};
    opts.nzkBL1 = man.nzkBL1.String{man.nzkBL1.Value};
    if strcmp(opts.nzkBL1,'Yes')
        opts.nzkBL1 = true;
    end
    opts.nzkSplCor = man.nzkSplCor.String{man.nzkSplCor.Value};
    if strcmp(opts.nzkSplCor,'Yes')
        opts.nzkSplCor = true;
    end
    
    opts.mzRes = [];
    
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

