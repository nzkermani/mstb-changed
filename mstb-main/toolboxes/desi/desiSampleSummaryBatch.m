function [numAnnos] = desiSampleSummaryBatch(fold,list)
%desiSampleSummaryBatch - run this for a bunch of MAT files

% Get a list of mat files in the folder
if nargin == 1
    allF = fileFinderAll(fold,'mat',true);
else
    allF = repmat({fold},[size(list,1) 2]);
    allF(:,2) = list;
end
numF = size(allF,1);

% Create folder one level further up...
[pf,sf] = previousFolder(fold);
saveF = [pf sf '-SummaryImages' filesep];
if ~exist(saveF,'dir')
    mkdir(saveF);
end

% Numbers of annotations...
%classes = {'Tumour';'Stroma';'Myometrium';'Background'};
%numAnnos = cell(numF,numel(classes)+1);
numAnnos = [];
mzVals = [726.544 739.5114 750.544 833.518 834.520 834.5315 887.569];

% Loop through each file - not all may be relevant!
for n = 1:numF
    
    % Prepare file name for saving...
    tmpF = allF{n,2}(1:end-4);
    newF = [saveF tmpF];
    if exist([newF '.png'],'file')
        continue;
    end
    
    numAnnos{n,1} = allF{n,2};

    % Open the file...
    if nargin == 1
        tmp = open([allF{n,1} filesep allF{n,2}]);
    else
        tmp = open([allF{n,1} allF{n,2} '.mat']);
    end
    
    % Does it have the correct bit?
    if ~isfield(tmp,'dpn')
        continue;
    end
    if ~isfield(tmp.dpn,'anno')
        continue;
    end
    [tmp.dpn.anno] = xxxAnnotationConsolidate(tmp.dpn.anno);
    
%     % Let's determine the number of annotations for each tissue type
%     [mask,histID,~] = desiAnnotationExtract(tmp.dpn);
%     fx = mask > 0;
%     these = histID(fx)
%     [un,~,ui] = unique(these)
%     for r = 1:numel(un)
%         
%         for s = 1:numel(classes)
%             tmpG = strfind(lower(un{r}),lower(classes{s}))
%             
%             if ~isempty(tmpG)
%                 numAnnos{n,s+1} = sum(ui == r);
%             end
%             
%         end
%         
%     end
    
    
    
    % Run the function
    try
        ff = desiSampleSummaryV2(tmp.dpn,mzVals);    
    
        % Save the figure
        pubFig(newF,'png');
    
        % Clsoe the figure
        close(ff.fig);
        
    catch err
        disp(allF{n,2});
        err
    end
    

end

%classes2 = [{'File'}; classes]'
%numAnnos = cat(1,classes2,numAnnos)

end

