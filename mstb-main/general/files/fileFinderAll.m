function [hits] = fileFinderAll(locn,type,omitFirst)
% fileFinderAll - determine all of the files of a specific type in a folder
% and subfolders

if strcmp(locn(end),filesep)
    locn = locn(1:end-1);
end

% Can we do it purely iteratively and dynamically, rather than going
% through all folders at the beginning...
[hits] = findFiles(locn,type,{'Path','Name'});

% Remove first header row, but carefully as other functions remove the
% first row manually...
if nargin == 2
    omitFirst = false;
else
    if ~islogical(omitFirst)
        omitFirst = false;
    end
end

if omitFirst
    hits = hits(2:end,:);
end
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hits] = findFiles(fold,type,hits)
% Find files of type in the current folder

% Start by generating the files within the 'root' folder
[root] = genFiles(fold);
numF = size(root,1);

%numH = size(hits,1);

for n = 1:numF
    
    % Check if it is a folder or not
    if root(n).isdir && ~strcmp(root(n).name(1),'.')
        
        % Then this is a folder. But don't look inside it if it is a raw
        % Waters file...
        if length(root(n).name) >= 4
            ext = fileExtension(root(n).name);
        else
            ext = '';
        end

        % Decide what to do
        if strcmp(ext,'raw') && strcmp(ext,type)            
            % Just take this raw file
            hits = cat(1,hits,{fold,root(n).name});
        else
            % Continue looking in this folder            
            [hits] = findFiles([fold filesep root(n).name],type,hits);
        end
        
    elseif length(root(n).name) > length(type)+1 ...
            && ~strcmp(root(n).name(1),'.') ...
            && ~strcmp(root(n).name(1),'~')
        
        % Then not a folder, so we ask if it is a file of correct type?
        %         idx = strfind(lower(root(n).name),['.' lower(type)]);
        %         if ~isempty(idx)
        %             hits = cat(1,hits,{fold,root(n).name});
        %         end          

        % Get the file extension
        ext = fileExtension(root(n).name);
        if strcmpi(type,ext)
            hits = cat(1,hits,{fold,root(n).name});
        end
        
            
        
    end
end
        
%hits = [];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [all] = genFiles(fold)
% Generate a list of files in 'fold'

all = dir(fold);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
