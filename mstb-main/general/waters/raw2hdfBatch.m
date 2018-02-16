function raw2hdfBatch
%raw2hdfBatch - process multiple files when I am away...

% Location of raw files
%rawLocn = 'E:\Box Sync\CRUK GC\Colorectal Xevo\Negative\';
%rawLocn = 'C:\Users\jmckenzi\Desktop\';
%rawLocn = 'Z:\Data\Jocelyn Waters Data\Raw data for Raman paper\Mouse brain\';
%rawLocn = 'Z:\Data\Xevo DESI\Colorectal Anna\';
rawLocn = 'Z:\Data\Xevo DESI\DESI ovarian 2ndBatch\raw\raw_neg\non_corrected\';

% Location in which to save the H5 files afterwards
%h5Locn = 'E:\Box Sync\H5 Convert\';
%h5Locn = 'C:\Users\jmckenzi\Desktop\';
%h5Locn = 'E:\Box Sync\H5 Convert\Pos\';
h5Locn = 'E:\Box Sync\H5 Convert\Ovarian-Neg\';

% Find all the raw files in the right place
allF = fileFinderAll(rawLocn,'raw');
allF = allF(6:end,:);
numF = size(allF,1);

% Remove all 'corr' from files and 'neg' too
%f1 = ~cellfun(@isempty,strfind(allF(:,2),'CORR'));
%f2 = ~cellfun(@isempty,strfind(allF(:,2),'corr'));
%f3 = ~cellfun(@isempty,strfind(allF(:,2),'neg'));
%fx = f1 | f2 | f3;
%fx = ~fx;
%allF = allF(fx,:);
%excl = [2 20 21 38 53 54 55 60];
%allF(excl,:) = [];
%numF = size(allF,1);

% Somewhere to store file names / errors
newFiles = cell(numF,2);

% Loop through
for n = 1:numF
    disp(allF{n,2});
    disp(int2str(n));
    [a,b] = raw2hdf(allF{n,1},allF{n,2},h5Locn);
    
    newFiles{n,1} = a;
    newFiles{n,2} = b;
    
end

assignin('base','newFiles',newFiles);

end

