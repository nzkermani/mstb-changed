% Dear James I wrote this code to test - stability (repeatability) of
% making intervals
% it has to main functions - the first one makeIntervals(path, file, ppm, mzRange)
%                          - the second testCoverage( path, file, ppm,
%                          intervalarray{i} ) 
% makeIntervals is the function that makes intervals
% testCoverage report the number of nonzeors mz in the original image and in the interval
% in addition the same info for intensities
% variable numberOfRuns gives you the option rerun makeIntervals as you
% wish. 

javaclasspath('packages/imzMLConverter/imzMLConverter.jar')
file = 'A62 CT S3-centroid.imzML';
path = '/Volumes/JSM/';
mzRange = [600 900];
ppm = 2.5;
interval = makeIntervalsv2(path, file, ppm, mzRange);
% test if the number of intervals is different among runs 
numberOfRuns = 10;
numbrOfIntervals = zeros(1,10);
% this is trick that I read somewhere to initiate a cell array.
intervalarray{numberOfRuns} = interval;
for i = 1:numberOfRuns
   intervalarray{i} = makeIntervalsv2(path, file, ppm, mzRange);
    [ mzInterval, spInterval ] = testCoverage( path, file, ppm, intervalarray{i} ,mzRange);
end
