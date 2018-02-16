function rawH5extractBatch(folder,files,list)
%rawH5extractBatch - gather the data from the H5 files according to the
%peak list determined in list (from rawH5global.m function)

% Get all files
%files = fileFinderAll(folder,'h5');
%files(1,:) = [];
numF = size(files,1);

% Generate a name for the h5 file
if ismac
    h5Name = ['/Users/jmckenzi/Documents/Box Sync/H5 Proc/Mac-H5-Test-PEAK-THRESH-2-' datestr(now,'yymmdd-HHMMSS') '.h5'];
else
    h5Name = ['E:\Box Sync\H5 Proc\Jocelyn' datestr(now,'yymmdd-HHMMSS') '.h5'];
end

% First let's save the mz vector that is common to all of the files
h5create(h5Name,'/mz',size(list.pks.mz));
h5write(h5Name,'/mz',list.pks.mz);

% Waitbar
wb = waitbar(0,'File Extraction');

for n = 1:numF
    
    % Create temporary pks structure
    pks.mz    = list.origMZ(:,n);
    pks.span  = list.pks.span';
    
    % Name of the file
    fName = [folder files{n,1}];
    
    % Extract function
    [data,mzImg] = rawH5extract(fName,pks);
    %data = rand(35,35,numel(pks.mz));
    
    % Need to save the 'data' and other things to an H5 file
    gpname = ['/file' int2str(n)];
    
    % Which bits are we going to save?
    
    h5create(h5Name,[gpname '/img'],size(data),...
        'Deflate',3,...
        'ChunkSize',[10 10 10]);
    h5write(h5Name,[gpname '/img'],data);
    
    h5create(h5Name,[gpname '/mzImg'],size(data),...
        'Deflate',3,...
        'ChunkSize',[10 10 10]);
    h5write(h5Name,[gpname '/mzImg'],mzImg);
    
    tic = nansum(data,3);
    h5create(h5Name,[gpname '/tic'],size(tic));
    h5write(h5Name,[gpname '/tic'],tic);
    
    h5writeatt(h5Name,gpname,'filePath',folder);
    h5writeatt(h5Name,gpname,'fileName',files{n,1});
        
    h5create(h5Name,[gpname '/rawMZ'],[size(data,3) 1]);
    h5write(h5Name,[gpname '/rawMZ'],list.origMZ(:,n));
    
    waitbar(n/numF,wb);   

end

delete(wb);


end

