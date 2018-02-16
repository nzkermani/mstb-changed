function [opts,flag] = desiGetProcOptions(file,force)
% desiGetProcOptions - determine the processing parameters (and defaults)
% for each of the various file types...

% Decide whether we can continue; if false then we should abort
flag = true;

% Options for processing imzML files 
opts.imzml.mzFrac   = 0.01;
opts.imzml.mzRes    = 0.002;
opts.imzml.ppmRes   = 8;
opts.imzml.mzRange  = [100 1000];
opts.imzml.nazanin  = 0;
opts.imzml.prompt = {'Fraction of empty pixels per m/z',...
    'm/z resolution',...
    'Maximum ±ppm shift',...
    'm/z range',...
    'Nazanin? (1 or 0)'};

% Options for processing RAW files
opts.raw.mzFrac     = 0.01;
opts.raw.mzRes      = 0.01;
opts.raw.ppmRes     = 20;
opts.raw.mzRange    = [100 1000];
opts.raw.peakDetect = 0;
opts.raw.recal      = 1;
opts.raw.prompt = {'Fraction of empty pixels per m/z',...
    'm/z resolution',...
    'Perform peak detection?',...
    'Maximum ±ppm shift',...
    'm/z range',...
    'Modify calibration function?'};

% Immediate return if 2 args in, as we just want the defaults and nothing
% else to be performed
if nargin == 2
    return
end
    

% Get the defaults in string-based format
xt = lower(file.ext);
switch xt
    case 'imzml'
        defs = {num2str(opts.(xt).mzFrac),...
            num2str(opts.(xt).mzRes),...
            int2str(opts.(xt).ppmRes),...
            num2str(opts.(xt).mzRange),...
            num2str(opts.(xt).nazanin)};
        
    case 'raw'
        defs = {num2str(opts.(xt).mzFrac),...
            num2str(opts.(xt).mzRes),...
            int2str(opts.(xt).peakDetect),...
            int2str(opts.(xt).ppmRes),...
            num2str(opts.(xt).mzRange),...
            int2str(opts.(xt).recal)};
        
    otherwise
        % Just return from here without fuss
        opts = [];        
        return
end

% Ask the user the question
resp = inputdlg(opts.(xt).prompt,'Processing parameters',1,defs);

if isempty(resp)
    % User aborted processing
    disp('User aborted the processing');
    opts = [];
    flag = false;
    return
end
    
% User defined options...
switch lower(file.ext)
    case 'imzml'
       
        % Convert to numbers - no error checking/validation though
        opts.(xt).mzFrac    = str2double(resp{1});
        opts.(xt).mzRes     = str2double(resp{2});
        opts.(xt).ppmRes    = str2double(resp{3});
        opts.(xt).mzRange   = str2num(resp{4}); %#ok<*ST2NM>
        opts.(xt).nazanin   = str2double(resp{5});

    case 'raw'

        % Convert to numbers - no error checking/validation though
        opts.raw.mzFrac     = str2double(resp{1});
        opts.raw.mzRes      = str2double(resp{2});
        opts.raw.peakDetect = str2double(resp{3});
        opts.raw.ppmRes     = str2double(resp{4});      
        opts.raw.mzRange    = str2num(resp{5});
        opts.raw.recal      = str2double(resp{6});
                
    otherwise
        % No need to take action
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%