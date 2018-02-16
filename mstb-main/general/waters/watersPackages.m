%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p1,p2] = watersPackages(filename)

libname = 'MassLynxRaw';

try
    if ~libisloaded(libname)
        loadlibrary(libname,'CMassLynxRawReader.h',...
            'addheader','MassLynxRawDefs');
    end
catch err
    error('Mex compilation issue / not Windows 64');
end

% Pointers to the file
p1 = calllib(libname,'newCMassLynxRawReader',filename);
p2 = calllib('MassLynxRaw','newCMassLynxRawScanReader',p1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
