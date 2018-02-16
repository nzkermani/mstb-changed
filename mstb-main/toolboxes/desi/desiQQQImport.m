function [masses_raw,X_raw,daughters_raw] = desiQQQImport(folder)
% desiQQQImport - read in QQQ data using Paolo's methodology. Note that
% this will only work on Windows machines with Visual Studio 2015
% installed. This will need to clarified and tested on various set ups.

% Simple Mac check
if ~ispc
    error('Cannot proceed on non-PC systems');
end

% Create a waitbar
wb = waitbar(0,'Importing QQQ');

% Path for the executable file - do not change this
exePath = [pwd filesep 'private\qqq-read_waters_raw.exe '];

% Be sure to remove the ending filesep if present
if strcmp(folder(end),filesep)
    folder = folder(1:end-1);
end

% Generate raw file list
list_raw = dir([folder, '\*.raw']);
list_raw = {list_raw.name};

% Read the raw data
num_samples = length(list_raw);

% Create a blank directory for temporary storage of data files
mkdir('./tmp');
masses = cell(num_samples, 1);
intensities = cell(num_samples, 1);
daughters = cell(num_samples, 1);
nscans = cell(num_samples, 1);
for i = 1:num_samples
    
    % delete ./tmp/*.*
    
    %fprintf('%d/%d\n', i, num_samples);
    
    raw_dir = ['"', folder, '\', list_raw{i}, '"'];
    
    system([exePath, raw_dir, ' ./tmp']);
    
    % Read the number of scans
    fid = fopen('./tmp/qqq-nscans.dat', 'r');
    nscans{i} = fread(fid, 'int');
    fclose(fid);
    
    fid = fopen('./tmp/qqq-intensities.dat', 'r');
    intensities{i} = fread(fid, 'float');
    fclose(fid);
    
    fid = fopen('./tmp/qqq-masses.dat', 'r');
    masses{i} = fread(fid, 'float');
    fclose(fid);
    
    fid = fopen('./tmp/qqq-daughters.dat', 'r');
    daughters{i} = fread(fid, 'float');
    fclose(fid);
    
    waitbar(i/num_samples,wb,['QQQ ' int2str(i) '/' int2str(num_samples)]);
    
end

clear i fid raw_dir

% Remove the temporary files
delete ./tmp/*.*
rmdir ./tmp

% Check that all the scans have the same number of masses
num_funcs_per_scan = cellfun(@(x) length(x), nscans);
if length(unique(num_funcs_per_scan)) ~= 1
    error('some scans contain different number of functions.');
end

% Transform the raw data into matrices (row x col x mass_number). Since the
% number of scans can vary from file to file, the final matrix will
% retrieve the minimum number of scans found across the files. All the
% additional scans are dropped. In this way the intensities can be easily 
% visualised.
num_funcs = unique(num_funcs_per_scan);
n_scans = min(cell2mat(nscans));
X_raw = nan(num_samples, n_scans, num_funcs);
masses_raw = nan(num_samples, n_scans, num_funcs);
daughters_raw = nan(num_samples, n_scans, num_funcs);
for i = 1:num_samples
    
    int_tmp = intensities{i};
    masses_tmp = masses{i};
    daughters_tmp = daughters{i};
    
    offset = 0;
    for j = 1:num_funcs
        
        X_raw(i, :, j) = int_tmp(offset+1:offset+n_scans);
        masses_raw(i, :, j) = masses_tmp(offset+1:offset+n_scans);
        daughters_raw(i, :, j) = daughters_tmp(offset+1:offset+n_scans);
        
        offset = offset + nscans{i}(j);
    end
    
    waitbar(i/num_samples,wb,'Quick');
    
end

%clear i j offset int_tmp masses_tmp daughters_tmp
delete(wb);
end