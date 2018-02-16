function bin = PPMPeakAlignmentv10parB1(mzTest,SpTest, ppm, binOfLengthOneOK,splitCorrection )
%% PPMPeakAlignmentv10par aligns m/z ratios using hierarchical agglomerative clustering 
% and ppm based heuristic rules for peak alignement  
%%
%    Input:  mzTest - [rows x columns x mz]
%            SpTest -  [rows x columns x intensity] 
%            ppm - part per million error of the instrument
%            binOfLengthOneOK - whether bin incluing only one ion is okay
%            (=1) or not (=0)
%    Output: bin - a structure that for each element of it includes information 
%                  of a group ions that represents the same ion species:
%                  mz - mz ratios 
%                  Spectra - intensities 
%                  centroidmz - centroid mzs
%                  centroidSpectra - centroid intensities 
%% Author: Nazanin z. Kermani, Imperial College London 2016.
% Reference: 1. Kazmi, S. A.,et al.Alignment of high resolution mass spectra: 
% Development of a heuristic approach for metabolomics. Metabolomics 2, 75–83 (2006).
%            2. to be published 


        %% Build super spectrum (/reference/target)
        % Error window is 2 times ppm 
        Eppm = 2*ppm/1e6;
        % pool all the spectra together and sort it out 
        [sorted_mz_value, sorted_mz_indx] = sort(mzTest(:));

        % spectra padded with lots of zero's, take the zero's out
        indx_zero = min(find(sorted_mz_value>0));

        % To facilitate paralle alignment 
        % caculate the distance between consecutive m/z values
        diff_mz = diff(sorted_mz_value(indx_zero:end));

        % Scale distances into 2 times ppm error window
        eppm = sorted_mz_value(indx_zero:end)*Eppm*2;

        % To facilitate paralle alignment 
        % find gaps in the spectra bigger that 2 times ppm error 
        % window (to create intervals)
        flag_indx = find(diff_mz > eppm(1:end-1));

        % creats intervals 
        number_of_intervals = size(flag_indx,1);
        count = 1;
        flag_indx = [0; flag_indx];

        %% Initiate intervals
        % Intervals partition the super spectrum and they can be aligned
        % independently. Interval is an array of structures that hold the
        % information of: 
        %                mz - m/z rations (ions) 
        %                I,J,Z - the coordinates on the 3d image cube 
        %                       (there is scope to delete this (TODO), not important but saves space) 
        %                spectra - the location on the initial sepectra 
        i=0;
        if(number_of_intervals>0)
            for i=1:(number_of_intervals)
                    interval(i).ions = sorted_mz_value((indx_zero+flag_indx(i)):(flag_indx(i+1)+indx_zero-1));
                    [interval(i).I, interval(i).J,  interval(i).Z] = ind2sub(size(mzTest)...
                        , sorted_mz_indx((indx_zero+flag_indx(i)):(flag_indx(i+1)+indx_zero-1))); 
                    interval(i).spectra = sorted_mz_indx((indx_zero+flag_indx(i)):(flag_indx(i+1)+indx_zero-1));
            end
         end
        interval(i+1).ions = sorted_mz_value((indx_zero+flag_indx(end)):end);
        [interval(i+1).I, interval(i+1).J,  interval(i+1).Z] = ind2sub(size(mzTest)...
            , sorted_mz_indx((indx_zero+flag_indx(end)):end));
        interval(i+1).spectra = sorted_mz_indx((indx_zero+flag_indx(end)):end);
        %% Do the alignment here in parallel for loop

        % intiate a parallel pool
        %numOfCPUs = getenv('NUMBER_OF_PROCESSORS');
        %poolObj = parpool('local', str2double(numOfCPUs));
%         poolObj = parpool('local');

       parfor i=1:(number_of_intervals+1)
        %for i=1:(number_of_intervals+1)
            % align the intervals individualy
            bin{i} = alignIntervalB1(interval(i),mzTest,SpTest, ppm, binOfLengthOneOK,splitCorrection) ;
        end
        % delete pool object
%         delete(poolObj);
        %% replace bin structure with bin of type cell to ease the downstream processing
        bin = binStruct2bin(bin);

        % clear all the variables except bin
        clearvars -except bin
end