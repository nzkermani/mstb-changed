function interval = makeIntervals(path, file, ppm, mzRange )
javaclasspath('packages/imzMLConverter/imzMLConverter.jar')
fullfile = [path file '.imzML'];
imzML = imzMLConverter.ImzMLHandler.parseimzML(fullfile);

% Get the image spectrum - note that this is perhaps the wrong function,
% and could use instead the more traditional one. I think that this
% requires options to be specified as well, which they are not here.
% Nona What options?
[mz,sp] = getMSImageCentroidMatrix(imzML);
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
    sorted_mz_value sorted_mz_indx TempSp dims
number_of_intervals = size(flag_indx,1);
flag_indx = [0; flag_indx];
%% Initiate intervals
% Intervals partition the super spectrum and they can be aligned
% independently. Interval is an array of structures that hold the
% information of:
%                mz - m/z rations (ions)
%                I,J,Z - the coordinates on the 3d image cube
%                       (there is scope to delete this (TODO), not important but saves space)
%                spectra - the location on the initial sepectra

% i=0, For the cases that there is no interval
i=0;

if(number_of_intervals>0)
    for i=1:(number_of_intervals)
        interval(i).ions = sorted_mz_value((flag_indx(i)+1):(flag_indx(i+1)));
        interval(i).spectra = sorted_mz_indx((flag_indx(i)+1):(flag_indx(i+1)));
        interval(i).intensity = TempSp(((flag_indx(i)+1):(flag_indx(i+1))));
    end
end
interval(i+1).ions = sorted_mz_value((flag_indx(end)+1):end);
interval(i+1).spectra = sorted_mz_indx((flag_indx(end)+1):end);
interval(i+1).intensity = TempSp(((flag_indx(end)+1):end));
disp([int2str(number_of_intervals) ' ---- Number of intervals for alignment']);
clearvars -except number_of_intervals interval ppm binOfLengthOneOK splitCorrection...
    method sorted_mz_indx sorted_mz_value sorted_mz_indx TempSp dims

end