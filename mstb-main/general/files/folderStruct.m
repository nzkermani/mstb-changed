function [ output_args ] = folderStruct( input_args )
%folderStruct - make folder structure from one place in another...

% Where the data is currently
hostA = '/Volumes/Untitled/Data Ottmar Golf/raw data/colorectal/DESI/CRC Set/raw data/';
hostB = {'DESI images 2012';'DESI images 2013';'DESI images 2014'};

% To where the data should be moved
target = '/Users/jmckenzi/Desktop/ML-ICL-Colo/';

% Which kinds of files to copy?
extn = {'imzML','ibd','ndpi'};
pics = {'png','tif','jpeg','jpg','xls','xlsx'};
if ~exist([target 'Misc/'],'dir')
    mkdir([target 'Misc/']);
end


% Loop through the folders creating, maybe copying the folders across...
numF = numel(hostB);
for n = 1:numF
    
    % Run the function...
    doCopy([hostA hostB{n}],target,extn,pics);
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doCopy(host,target,extn,pics)
% This is the function that does the donkey work

% Generate the folder structure in here... (not recursive)
all = dir(host);

numF = size(all,1);


for n = 1:numF
    
    % Check that it is suitable
    if all(n).isdir && length(all(n).name) > 3 && strcmpi(all(n).name(1),'a')
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
            copyfile(fullPath,newFold);            
        end
        
        % Find all image files and other perhaps useful things
        cmp2 = strcmp(pics,chk);
        if sum(cmp2) > 0 && ...
                ~strcmp(cont(r).name(1:4),'Data') && ...
                ~strcmp(cont(r).name(1:4),'MSIm') 
            
            % Full path of original and then copy (as above)
            fullPath = [host filesep all(n).name filesep cont(r).name];
            %disp(cont(r).name);
            copyfile(fullPath,[target 'Misc/']);
        end
    end
end







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
