function [info] = mtspImportTxtDump
%mtspImportTxtDump - get a list of the annotations for each of the files in
% the main folder

% Main directory?
fold = '/Volumes/JSM/DB/Metaspace/NewEngineDump/Raw/';

% Get a list of folders in this folder
allF = folderList(fold);
numF = size(allF,1);

% For each of these, get the file lists...
info = cell(numF,5);
for n = 1:numF
    
    % Just the file name
    info{n,1} = allF{n,1};
    
    % Get the names/paths of the text files in each folder
    subF = fileFinderAll([fold info{n,1}],'txt',true);
    info{n,2} = subF;
    
    % Now we need to split some of this up
    [info{n,3},info{n,4},info{n,5}] = extractInfo(subF(:,2));
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [idx,anno,adct] = extractInfo(names)
% Extract info from the file names

und = strfind(names,'_');

sz = size(names,1);

idx = zeros(sz,1);
anno = cell(sz,1);
adct = cell(sz,1);

for n = 1:sz
    
    ii = und{n,1};
    
    idx(n,1) = str2double(names{n,1}(ii(1)+1:ii(2)-1));
    anno{n,1} = names{n,1}(ii(2)+1:ii(3)-1);
    adct{n,1} = names{n,1}(ii(3)+1:end-4);
    
end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%