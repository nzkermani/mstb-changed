function [ output_args ] = dbChangeHistID
%dbChangeHistID - nuclear option for permanently modifying the histID of a
%series of files.  Not for the faint hearted.

% Supply the path here, which contains the mat files to be changed
path = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Final/Lymph/';

% Supply the translation file here
trans = '/Volumes/JSM/DB/LymphRecalNonMTSP/Lymph Histology/LymphTranslation.xlsx';

% Simple checks...
if validateFcn(path,trans)
    return
end

% First we need to get the translations
[~,~,c] = xlsread(trans);
c = c(:,1:2);

% Get a list of all mat files in the path
allF = fileFinderAll(path,'mat',true);
numF = size(allF,1);

% Loop through each file
for n = 1:numF
    
    % Determine save name
    saveName = [allF{n,1} '-Translated' filesep allF{n,2}];
    if exist(saveName,'file')
        disp([int2str(n) '---SKIP']);
        continue;
    end
    
    % Open the file
    tmp = open([allF{n,1} filesep allF{n,2}]);
    if ~isfield(tmp,'dpn')
        continue;
    end
    
    % Get the histID
    oh = tmp.dpn.anno(:,5);
    nh = oh;
    
    % Unique...
    [unq,~,ind] = unique(oh);
    
    % Translate
    for r = 1:numel(unq)
        
        % Compare to c list
        cmp = strcmp(c(:,1),unq{r});
        
        % Replace
        fx = ind == r;
        nh(fx,:) = c(cmp,2);
        
    end
    
    % Now with the translations, replace in oanno
    tmp.dpn.anno(:,5) = nh;
    
    % Consider saving the files - let's test first
    dpn = tmp.dpn;
    
    save(saveName,'dpn');
        
    disp([int2str(n) '/' int2str(numF) ' --- ' allF{n,2}]);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [check] = validateFcn(path,trans)
% Check that the path/file exist

if exist(path,'dir') && ...
        exist(trans,'file') && ...
        strcmp(getUser,'jmckenzi')    
    check = false;
else
    check = true;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%