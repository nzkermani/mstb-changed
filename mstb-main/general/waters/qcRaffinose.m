function [ imgs ] = qcRaffinose(fold)
%qcRaffinose - extrac images for lock mass corrected files showing the
%intensities of raffinose

% Raf m/z
mz = 503.1627;
ppm = 0.5;

% Get the files
allF = fileFinderAll(fold,'raw',true);

% Strip out the ones we don't want, i.e. only want corr503 in neg
%f1 = ~cellfun(@isempty,strfind(allF(:,2),'NEG'));
%f2 = ~cellfun(@isempty,strfind(allF(:,2),'neg'));

%f3 = ~cellfun(@isempty,strfind(allF(:,2),'CORR'));
%f4 = ~cellfun(@isempty,strfind(allF(:,2),'corr'));

%fx = (f1 | f2);
%fy = (f3 | f4);

%fz = fx & fy;

%allF = allF(fz,:);

%return

numF = size(allF,1);

% Loop through the files
imgs = cell(numF,2);
for n = 1:numF
    
    disp(int2str(n));
    
    try
        
    % Is the file already copied to this computer?
%     newFile = ['E:\Data\Olivia\' allF{n,2}];
%     if ~exist(newFile,'file')  
%         copyfile([allF{n,1} filesep allF{n,2}],newFile);
%     end
     
    % Import the data...
    [sp,numP,xy,xy2D] = desiReadRaw([allF{n,1} filesep allF{n,2}],true);
    
    % Create the image
    imgs{n,1} = rawImageMake(sp,xy2D,mz,ppm);
    imgs{n,2} = allF{n,2};
    
    catch
        imgs{n,1} = rand(100,100);
        imgs{n,2} = allF{n,2};
    end
    
end

end

