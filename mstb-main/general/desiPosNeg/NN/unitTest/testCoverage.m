function [ mzInterval, spInterval ] = testCoverage( path, file,ppm, interval, mzRange );
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
javaclasspath('packages/imzMLConverter/imzMLConverter.jar')
fullfile = [path file];
imzML = imzMLConverter.ImzMLHandler.parseimzML(fullfile);

% Get the image spectrum - note that this is perhaps the wrong function,
% and could use instead the more traditional one. I think that this
% requires options to be specified as well, which they are not here.
% Nona What options?
[mz,sp] = getMSImageCentroidMatrix(imzML,mzRange);
%[mz,sp] =getMSImageCentroidMatrixRange(imzML,[600 1000]);
clear imzML
% here can go the data transformation(log, sqrt,...) and 
% normalization(TIC,...) if necessary
dims = size(sp);

%% Build super spectrum (/reference/target)
% Error window is 2 times ppm
Eppm = 2*ppm/1e6;
% pool all the spectra together and sort it out
[sorted_mz_value, sorted_mz_indx] = sort(mz(:));
clear mz
sorted_mz_value = sparse(sorted_mz_value);
indx_zero = min(find(sorted_mz_value>0));
sorted_mz_value(1:indx_zero) = [];
sorted_mz_indx(1:indx_zero) = [];
TempSp = sp(:);
clear sp
TempSp = TempSp(sorted_mz_indx);
% spectra padded with lots of zero's, take the zero's out
% To facilitate paralle alignment
% caculate the distance between consecutive m/z values
diff_mz = diff(sorted_mz_value);

% Scale distances into 2 times ppm error window
eppm = sorted_mz_value*Eppm*2;

% To facilitate paralle alignment
% find gaps in the spectra bigger that 2 times ppm error
% window (to create intervals)
flag_indx = find(diff_mz > eppm(1:end-1));

% creats intervals
clearvars -except TempSp sorted_mz_value  sorted_mz_indx flag_indx ppm...
    binOfLengthOneOK splitCorrection method ...
    sorted_mz_value sorted_mz_indx TempSp dims interval
number_of_intervals = size(flag_indx,1);
flag_indx = [0; flag_indx];

%% reverse making intervals
i=0;
mzInterval = zeros(1,length(TempSp));
spInterval = zeros(1,length(TempSp));
if(number_of_intervals>0)
    for i=1:(number_of_intervals)
        mzInterval((flag_indx(i)+1):(flag_indx(i+1))) = interval(i).ions;
        spInterval((flag_indx(i)+1):(flag_indx(i+1))) =  interval(i).intensity; 
    end
end
mzInterval((flag_indx(end)+1):end) = interval(i+1).ions; 
spInterval((flag_indx(end)+1):end) = interval(i+1).intensity; 

disp([int2str(length(mzInterval)) ', ' int2str(length(sorted_mz_value)) '---- Number of mz stored in the intervals for alignment, Number of non zetro mz to begin with']);

disp([int2str(length(spInterval)) ', ' int2str(length(TempSp)) ' ---- Number of intesities stored in the intervals for alignment, Number of non zetro intesities to begin with']);

% now test if they are equal
if(all(spInterval' == TempSp))
    disp([int2str(all(spInterval' == TempSp)) ' ---- interval intesity is the same as the original image']);
else
    disp([int2str(all(spInterval' == TempSp)) ' ---- DANGER, interval intesity is not the same as the original image']);
end

% now test if they are equal
if(all(mzInterval == sorted_mz_value'))
    disp([int2str(all(mzInterval == sorted_mz_value')) ' ---- interval mz is the same as the original image']);
else
    disp([int2str(all(mzInterval == sorted_mz_value')) ' ---- DANGER, interval mz is not the same as the original image']);
end
    
clearvars -except mzInterval spInterval

end

