function [ allF ] = folderStruct2( input_args )
%folderStruct - moves files from one place according to their type

% Where the data is currently
hostA = '/Volumes/Untitled/Data Ottmar Golf/raw data/colorectal/DESI/CRC Set/raw data/';
hostB = {'DESI images 2012';'DESI images 2013';'DESI images 2014'};

hostA = '/Volumes/Data/Data/Liver/';
hostB = {'3D-Liver'};

% To where the data should be moved
target = '/Users/jmckenzi/Desktop/ML-ICL-3DLiver/';

% Which kinds of files to copy?
extn = {'imzML','ibd','ndpi','jpeg','jpg','png'};%,'tif','jpeg','jpg'}


% Loop through the folders creating, maybe copying the folders across...
numF = numel(hostB);
for n = 1:1%numF
    
    % Run the function...
    [allF] = doCopy([hostA hostB{n}],target,extn);
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [allF] = doCopy(host,target,extn)
% This is the function that does the donkey work

% Generate the folder structure in here... (not recursive)
all = dir(host);

numF = size(all,1);

i = 0;
allF = cell(200,3);

for n = 1:numF
    
    % Check that it is suitable
    if all(n).isdir && length(all(n).name) > 3 %&& strcmpi(all(n).name(1),'a')
        % Continue
    else
        % Break through this loop
        continue;
    end
    
    % First we copy the folder to the new place...
    newFold = [target all(n).name];
    if ~exist(newFold,'dir')
        mkdir(newFold);
    end
    
    % File list of contained files
    cont = dir([host filesep all(n).name]);
    
    % Now we could look to copy the important files...
    for r = 1:size(cont,1)
        
        % Get the extension
        chk = getExt(cont(r).name);
        
        % Compare to master extension list
        cmp = strcmp(extn,chk);
        
        if sum(cmp) > 0 && ...
                ~strcmp(cont(r).name(1:4),'test') && ...
                ~strcmp(cont(r).name(1:4),'MSIm') && ...
                ~strcmp(cont(r).name(1:4),'Data')
            
            % Full path of original
            fullPath = [host filesep all(n).name filesep cont(r).name];
            disp(cont(r).name);
            
            % Copy this file...
            %copyfile(fullPath,[target all(n).name]);
            
            % Megalist of all files
            i = i + 1;
            allF{i,1} = [host filesep all(n).name];
            allF{i,2} = cont(r).name;
            allF{i,3} = chk;
            
        end
        
    end
end

allF = allF(1:i,:);






end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ext] = getExt(name)
% Return the extension of the file

try
    dots = strfind(name,'.');
    ext = name(dots(end)+1:end);
catch
    ext = '';
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
