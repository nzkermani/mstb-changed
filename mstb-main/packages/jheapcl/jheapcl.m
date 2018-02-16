function [] = jheapcl(verbose)

if nargin < 1
    verbose = 0;
end

try
    org.dt.matlab.utilities.JavaMemoryCleaner.clear(verbose);
catch
    javaaddpath(which('MatlabGarbageCollector.jar'));
    org.dt.matlab.utilities.JavaMemoryCleaner.clear(verbose)
end

% Use this for silent cleanup
% org.dt.matlab.utilities.JavaMemoryCleaner.clear(1)

% Decomment this for verbose cleanup
% org.dt.matlab.utilities.JavaMemoryCleaner.clear(1)

end