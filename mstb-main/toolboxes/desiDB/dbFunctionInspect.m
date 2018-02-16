function dbFunctionInspect(~,~,data,opts,fig,pan)
%dbFunctionInspect - do the processing and manipulation of the data

% Let's strip out the deleted groups first
[data] = removeSpectra(data,opts);

% Now we need to combine spectra according to the options. Note that we
% overwrite the original data because of potential size issues
[data] = specCombine(data,opts);

% Here is where we add the new histological information to the data, based
% on the specified Excel spreadsheet. We don't want t remove things, just
% change the labels of the groups afterwards
[op] = dbAddHistID(data.meta.histID,pan);
if ~isempty(op)
    data.meta.newHist = op;
end

% Add in a separate part...
[unq,~,~] = unique(data.meta.histID);
fx = ~cellfun(@isempty,strfind(unq,'Myometrium_'));
if sum(fx) > 0
    newHist = repmat({'Other'},size(data.meta.histID));
    
    fx = ~cellfun(@isempty,strfind(data.meta.histID,'Myometrium'));
    newHist(fx) = {'Myometrium'};
    
    fx = ~cellfun(@isempty,strfind(data.meta.histID,'Stroma'));
    newHist(fx) = {'Stroma'};
    
    fx = ~cellfun(@isempty,strfind(data.meta.histID,'Tumour'));
    newHist(fx) = {'Tumour'};
    
    data.meta.newHist = newHist;
    
end
    

data.opts = opts;
data.path = fig.fig.Name;

assignin('base','msa',data);

disp('-----------------------------------------');
disp('Finished!');
disp('Now you need to run: stats');
disp('Then import the variable called msa');
disp('-----------------------------------------');

%dbStats(data);
%launchStats(data.mz,data.sp,data.meta);

disp('Done');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = removeSpectra(data,opts)
% Remove unwanted groups from the structure, and make the new vector for
% histID

% How many original groups?
numG = size(opts.names,1);

% New group information
newHist = cell(size(data.meta.histID));

% Find the old names for the newHist vector
for n = 1:numG
        
    fx = strcmp(data.meta.histID,opts.names{n,1});    
    
    if strcmp(opts.names{n,2},'Original')
        newHist(fx,1) = opts.names(n,1);
    else
        newHist(fx,1) = opts.names(n,2);
    end
end

% Now find the entries marked as 'Delete'
fx = strcmpi(newHist,'delete');

% And strip out the relevant bits
data.sp(fx,:) = [];
if isfield(data.meta,'date')
    data.meta.date(fx,:)   = [];
end
data.meta.histID(fx,:) = [];
data.meta.fileID(fx,:) = [];
data.meta.newHist = newHist(~fx,:);
if isfield(data.meta,'foldID')
    data.meta.foldID = data.meta.foldID(~fx,:);
end
if isfield(data.meta,'foldID2')
    data.meta.foldID2 = data.meta.foldID2(~fx,:);
    data.meta.foldID3 = data.meta.foldID3(~fx,:);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newData] = specCombine(data,opts)
% Combine the spectra by file and tissue type

% No need to combine if we want all pixels
if strcmpi(opts.pixels,'all')
   
    newData.mz = data.mz;
    newData.sp = data.sp;
    newData.meta.fileID = data.meta.fileID;
    newData.meta.histID = data.meta.newHist;
    
    if isfield(newData.meta,'date')
        newData.meta.date   = data.meta.date;
    end

    try
        newData.meta.foldID = data.meta.foldID;
    catch
    end
    
    try
        newData.meta.foldID2 = data.meta.foldID2;
        newData.meta.foldID3 = data.meta.foldID3;
    catch
    end
    
    % Metaspace specific things...
    if isfield(data,'annos')
        
        % Include the annotations
        newData.annos = data.annos;
        
        % Let's try to convert the annotations into m/z values, to make
        % life make more sense...
        tmpMZ = mtspForm2MZ(newData.annos);        
        [newData.mz,idx] = sort(tmpMZ);
        newData.sp = newData.sp(:,idx);
        newData.annos = newData.annos(idx,:);
        
        % Determine the quantity of annotations
        numAnno = round(100 * sum(newData.sp > 0,2) / size(newData.sp,2));
        numAnno = round(numAnno);
        newData.meta.numAnno = zeros(size(numAnno,1),1);
        for n = 1:size(numAnno,1)
            newData.meta.numAnno(n,1) = numAnno(n,1);%int2str(numAnno(n));
        end
        
    end    

    
    return
end

% Unique files and tissue types
[unqF,firstF,indF] = unique(data.meta.fileID);
numF = numel(unqF);

% Store the combined data here
newData = struct('sp',[],'histID',[],'fileID',[]);
i = 0;

% Loop through each file
for n = 1:numF
    
    % Indices for file n
    fx = indF == n;
    
    % Temporary spectra / histIDs
    tmpSp = data.sp(fx,:);
    tmpHi = data.meta.newHist(fx,:);
    
    % Date for this file?
    try
        fileDate = data.meta.date(firstF(n));
    catch
        fileDate = now;
    end
    
    % What is the subType / foldID?
    try
        foldID = data.meta.foldID{firstF(n)};
    catch
        foldID = 'Unknown';
    end
    try
        foldID2 = data.meta.foldID2{firstF(n)};
        foldID3 = data.meta.foldID3{firstF(n)};
    catch
        foldID2 = 'Unknown';
        foldID3 = 'Unknown';
    end
    
    % Find unique classes in tmpHi
    [unqH,~,indH] = unique(tmpHi);
    numH = numel(unqH);
    
    % Loop through each class
    for r = 1:numH
        
        % Indices
        fy = indH == r;
        
        % Increase the counter
        i = i + 1;
        
        % Take some spectra...
        switch lower(opts.pixels)            
            case {'tissue type mean','mean spectrum'}
                newData(i).sp = nanmean(tmpSp(fy,:),1);
                newData(i).histID = {unqH{r}};
                newData(i).fileID = {unqF{n}};
                newData(i).date   = fileDate;
                newData(i).foldID = {foldID};
                newData(i).foldID2 = {foldID2};
                newData(i).foldID3 = {foldID3};
                
            case {'tissue type median','median spectrum'}
                newData(i).sp = nanmedian(tmpSp(fy,:),1);
                newData(i).histID = {unqH{r}};
                newData(i).fileID = {unqF{n}};
                newData(i).date   = fileDate;
                newData(i).foldID = {foldID};
                newData(i).foldID2 = {foldID2};
                newData(i).foldID3 = {foldID3};

            case {'random','random 5'}
                numPix = min([5 sum(fy)]);
                randInd = randperm(sum(fy),numPix);
                
                newData(i).sp = tmpSp(randInd,:);
                newData(i).histID = repmat({unqH{r}},[numPix 1]);
                newData(i).fileID = repmat({unqF{n}},[numPix 1]);
                newData(i).date   = repmat(fileDate,[numPix 1]);
                newData(i).foldID = repmat({foldID},[numPix 1]);
                newData(i).foldID2 = repmat({foldID2},[numPix 1]);
                newData(i).foldID3 = repmat({foldID3},[numPix 1]);

            otherwise
                % Do nothing                
        end
        
    end
    
end

tmp = newData;
clear newData
newData.mz = data.mz;
newData.sp = vertcat(tmp.sp);
newData.meta.fileID = vertcat(tmp.fileID);
newData.meta.histID = vertcat(tmp.histID);
newData.meta.date   = vertcat(tmp.date);
newData.meta.foldID = vertcat(tmp.foldID);
newData.meta.foldID2 = vertcat(tmp.foldID2);
newData.meta.foldID3 = vertcat(tmp.foldID3);

if isfield(data,'annos')
    
    % Save the annotations...
    newData.annos = data.annos;

    % Determine the new m/z values
    tmpMZ = mtspForm2MZ(newData.annos);        
    [newData.mz,idx] = sort(tmpMZ);
    newData.sp = newData.sp(:,idx);
    newData.annos = newData.annos(idx,:);    
    
    % Determine the quantity of annotations
    numAnno = sum(newData.sp > 0,2);
    %numAnno = round(numAnno);
    newData.meta.numAnno = zeros(size(numAnno,1),1);
    for n = 1:size(numAnno,1)
        newData.meta.numAnno(n,1) = numAnno(n,1);%int2str(numAnno(n));
    end


end    

% Add these to the meta arguments...
%newData.meta.list1 = alt1;
%newData.meta.list2 = alt2;



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alt1,alt2] = additionalSplit(meta)

% Find another way to split up the metadata according to unique words
numM = size(meta.histID,1);
ff = [meta.histID repmat({' '},[numM 1])]';
ff = [ff{:}];

% Find spaces
spaces = [0 strfind(ff,' ')];
idx = [(spaces(1:end-1) + 1)' (spaces(2:end) - 1)'];
allW = cell(size(idx,1),1);
for n = 1:size(idx,1)
    allW{n,1} = ff(idx(n,1):idx(n,2));
end

% Separate according to alternative words
aw1 = unique(allW(1:2:end-1));
aw2 = unique(allW(2:2:end));

alt1 = repmat({'NA'},[numM 1]);
for n = 1:numel(aw1)    
    fx = ~cellfun(@isempty,strfind(meta.histID,aw1{n}));
    alt1(fx) = aw1(n);
end

alt2 = repmat({'NA'},[numM 1]);
for n = 1:numel(aw2)    
    fx = ~cellfun(@isempty,strfind(meta.histID,aw2{n}));
    alt2(fx) = aw2(n);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alt1,alt2] = oliviaSplit(meta)

classes = {'Background','Myometrium','Stroma','Tumour'};
numC = numel(classes);

numM = size(meta.histID,1);
newHist = repmat({'Other'},size(meta.histID));

for n = 1:numC
    
    fx = strfind(meta.histID,classes{n});
    fx = ~cellfun(@isempty,fx);
    
    newHist(fx,1) = classes(n);
    
end

alt1 = newHist;
alt2 = newHist;

return



% Find another way to split up the metadata according to unique words
ff = [meta.histID repmat({' '},[numM 1])]';
ff = [ff{:}];

% Find spaces
spaces = [0 strfind(ff,' ')];
idx = [(spaces(1:end-1) + 1)' (spaces(2:end) - 1)'];
allW = cell(size(idx,1),1);
for n = 1:size(idx,1)
    allW{n,1} = ff(idx(n,1):idx(n,2));
end

% Separate according to alternative words
aw1 = unique(allW(1:2:end-1));
aw2 = unique(allW(2:2:end));

alt1 = repmat({'NA'},[numM 1]);
for n = 1:numel(aw1)    
    fx = ~cellfun(@isempty,strfind(meta.histID,aw1{n}));
    alt1(fx) = aw1(n);
end

alt2 = repmat({'NA'},[numM 1]);
for n = 1:numel(aw2)    
    fx = ~cellfun(@isempty,strfind(meta.histID,aw2{n}));
    alt2(fx) = aw2(n);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
