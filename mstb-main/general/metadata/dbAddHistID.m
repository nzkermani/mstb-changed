function [ newHist ] = dbAddHistID(histID,fig)
%dbAddHistID - add new information to the data based on an excel
%spreadsheet of old and new histID

% Get the file name...
file = fig.translate.UserData;
if isempty(file)
    newHist = [];
    return
end

% First we need to get the translations
[~,~,c] = xlsread([file.dir file.nam]);
c = c(:,1:2);

% What are the unique histIDs from the old data?
[unqH,~,indH] = unique(histID);
numH = numel(unqH);

% Duplicate existing histID - use existing hists where no new hist is
% provided.
newHist = histID;

% Work through each of the existing histIDs and replace where necessary
for n = 1:numH
    
    % This histID
    th = unqH{n};
    fx = indH == n;
    
    % Can we find it in the translation list?
    trl = strcmp(c(:,1),th);
    
    % Only continue if we find its entry; if we don't then we just keep the
    % existing value...
    if sum(trl) ~= 1
        continue;
    end
    
    % If yes, what is its translation?
    trans = c(trl,2);
    
    % Check to see that it isn't numeric, i.e. 'NaN'
    if isnumeric(trans{1})
        % Do nothing...
        continue;
    end
    
    % Now just dump new version in newHist
    newHist(fx,1) = trans;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
