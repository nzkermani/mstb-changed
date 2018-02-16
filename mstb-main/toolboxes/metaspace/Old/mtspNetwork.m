function [op] = mtspNetwork(allA)
%mtspNetowrk - uses the output from the desiDB functions, which is
%typically called the 'msa' variable, i.e. the input to the stats toolbox.

% Where are the files? Do we want to use all of them or just some of them?
path = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Final/';
%path = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Final/Breast/';

% Have the annotations been provided?
if nargin == 0
    allA = [];
end
if isempty(allA)
    [allA] = getAnnos(path);
    op.allA = allA;
    return
else
    op.allA = allA;
end

% Convert annotations into one long cell
af = vertcat(allA.file);
aa = classMany2One(vertcat(allA.anno),'@');
at = vertcat(allA.tissue);

[op.gr,op.conn] = nwObs(af,aa,at);

% Determine the total number of annotations per observation
totAnno = zeros(size(op.conn,1),1);
for n = 1:size(totAnno,1)
    totAnno(n,1) = size(op.allA(n).anno,1);
end
op.gr.Nodes.NumAnno = totAnno;

% Now we need to link the histological annotations to the correct
% observations, noting that the orders are no longer the same, so we need
% to do this cleverly
allH = unique(vertcat(allA.histID));
for n = 1:numel(allH)
    
    % Empty matrix to be either 0 or 1 for absence/presence
    tmp = zeros(size(op.conn,1),1);

    % Loop through each observation
    for r = 1:size(op.allA,2)
        
        % What is this observation called?
        thisFile = op.allA(r).file{1};
        thisHist = op.allA(r).histID;
        
        % Find the index in op.gr.Nodes.Name
        comp = strcmp(op.gr.Nodes.Name,thisFile);
        
        % Compare to what we are actually looking for, i.e. allH{n}
        cmpHist = strcmp(thisHist,allH{n});
        
        % Add to tmp
        if any(cmpHist)
            tmp(comp,1) = 1;
        end
        
    end
    
    % Format name to be compatible
    nam = allH{n};
    vec = isstrprop(nam,'alphanum');
    nam(~vec) = '_';
    disp(nam);

    % Add as node label
    op.gr.Nodes.(nam) = tmp;
    
end             

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [allAnno] = getAnnos(path)
  
% Get the list of files here
allF = fileFinderAll(path,'mat',true);

% How many files are there?
numF = size(allF,1);

% Empty structure in which to store all annotations
allAnno = struct('file',[],'anno',[],'tissue',[],'histID',[]);

% Loop through the files reading just the annotations
for n = 1:numF
    
    % Read in
    tmp = open([allF{n,1} filesep allF{n,2}]);
    a = tmp.dpn.mtsp.annos;
    
    % What is the containing folder called?
    [~,contFold] = previousFolder(allF{n,1});
    
    % Save to structure
    allAnno(n).file = repmat(allF(n,2),[size(a,1) 1]);
    allAnno(n).anno = a;
    allAnno(n).tissue = repmat({contFold},[size(a,1) 1]);
    allAnno(n).histID = unique(tmp.dpn.anno(:,5));
    disp([int2str(n) '/' int2str(numF)]);
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gr,conn] = nwObs(file,anno,tissue)
% Link observations according to the annotations

% Find all unique annotations
[unqA,~,indA] = unique(anno);
numA = numel(unqA);

% Find all unique files
[unqF,firstF,~] = unique(file);
numF = numel(unqF);

% File connectivity matrix
conn = zeros(numF,numF);

%figure;

% Now run through each annotation and find links across the files
for n = 1:numA
    
    % Indices of this node
    fx = indA == n;
    
    % Which files have this annotation?
    j = unique(file(fx,1));
    
    if size(j,1) > 0
        % Ignore annotations in a single file
        
        % Find connections
        ism = ismember(unqF,j);
        
        % Add to the conn matrix (don't worry about the diagonal)
        conn(ism,ism) = conn(ism,ism) + 1;
        
        %imagesc(conn)
        %pause(0.1);
            
    end
        
end

% Remove diagonal elements
dg = eye(numF);
conn(dg == 1) = 0;

% Create the graph!
gr = graph(conn,unqF,'upper');

% Add various labels to the graph
gr.Nodes.Tissue = tissue(firstF);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gr] = fileGraph(all)
% Make a table showing files linked by annotations

% Find all unique annotations
[unqA,~,indA] = unique(all(:,3));
numA = numel(unqA);

% Find all unique files
[unqF,firstF,~] = unique(all(:,2));
numF = numel(unqF);

% File connectivity matrix
conn = zeros(numF,numF);

% Now run through each annotation and find links across the files
for n = 1:numA
    
    % Indices of this node
    fx = indA == n;
    
    % Which files have this annotation?
    j = unique(all(fx,2));
    
    if size(j,1) > 1
        % Ignore annotations in a single file
        
        % Find connections
        ism = ismember(unqF,j);
        
        % Add to the conn matrix (don't worry about the diagonal)
        conn(ism,ism) = conn(ism,ism) + 1;
    end
        
end

% Now determine a lot of the other interesting things
charge = all(firstF,end);
% if numel(chargeTypes) == 2    
%     charge = all(firstA,end);
%     [conn] = linkCharges(conn,chemForm,charge);
% end

% Remove diagonal elements
dg = eye(numF);
conn(dg == 1) = 0;

% Create the graph!
gr = graph(conn,unqF,'upper');

breast = ~cellfun(@isempty,strfind(unqF,'ICL//BRB'));
colo = ~cellfun(@isempty,strfind(unqF,'ICL//A'));
astra = ~cellfun(@isempty,strfind(unqF,'AstraZeneca'));
genen = ~cellfun(@isempty,strfind(unqF,'Genen'));
liver = ~cellfun(@isempty,strfind(unqF,'TopRight'));
oesLN = ~cellfun(@isempty,strfind(unqF,'ICL//LN'));
oesTu = ~cellfun(@isempty,strfind(unqF,'ICL//TO'));

tissue = repmat({'Ovarian',},[numel(unqF),1]);
tissue(breast) = {'Breast'};
tissue(colo) = {'Colorectal'};
tissue(liver) = {'Liver'};
tissue(astra) = {'Astra'};
tissue(genen) = {'Genentech'};
tissue(oesLN) = {'OES-Lymph'};
tissue(oesTu) = {'OES-Tumour'};

gr.Nodes.Tissue = tissue;
gr.Nodes.Charge = charge;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gr] = annoGraph(all)
% Make a table showing unique annotations

% Find all unique annotations
[unqA,firstA,indA] = unique(all(:,3));
numA = numel(unqA);
mzAnno = all(firstA,5);

% Find all unique files
[unqF,~,indF] = unique(all(:,2));
numF = numel(unqF);

% File connectivity matrix
conn = zeros(numA,numA);

% Now run through each annotation and find links across the files
for n = 1:numF
    
    % Indices of this node
    fx = indF == n;
    
    % Which files have this annotation?
    j = unique(all(fx,3));
    
    if size(j,1) > 1
        % Ignore annotations in a single file
        
        % Find connections
        ism = ismember(unqA,j);
        
        % Add to the conn matrix (don't worry about the diagonal)
        conn(ism,ism) = conn(ism,ism) + 1;
    end
        
end

% For each node, let's work down the list of chemical formulae and decide
% what name would be better than the simple chemical formula
chemForm = all(firstA,3);
lmDB = open('/general/lipid/DB1.mat');
[chemName] = mtspFormulaTranslate(chemForm,lmDB.lm);

% Format the chemical name to have the class instead...
isbr = cellfun(@strfind,chemName,repmat({'('},[size(chemName,1) 1]),...
    'UniformOutput',false);
isbr2 = find(~cellfun(@isempty,isbr));
chemClass = repmat({'Unknown'},[size(chemName,1) 1]);
for n = 1:numel(isbr2)    
    i = isbr2(n);    
    chemClass{i} = chemName{i}(1:isbr{i}-1);    
end

% Now determine a lot of the other interesting things
% chargeTypes = unique(all(:,end));
% if numel(chargeTypes) == 2    
%     charge = all(firstA,end);
%     [conn] = linkCharges(conn,chemForm,charge);
% end
    

% Remove diagonal elements
dg = eye(numA);
conn(dg == 1) = 0;

% Find values of 0 and just remove them - these annotations just get in the
% way
% fx = sum(conn,1) > 10;
% conn = conn(fx,fx);
% unqA = unqA(fx,:);
% firstA = firstA(fx,:);
% mzAnno = mzAnno(fx);
% chemName = chemName(fx);
% chemClass = chemClass(fx);
% charge = charge(fx);


% Create the graph!
gr = graph(conn,unqA,'upper');

% Add various characteristics to each node
gr.Nodes.mz         = mzAnno;
gr.Nodes.Chemical   = chemName;
gr.Nodes.Class      = chemClass;
%gr.Nodes.Charge     = charge;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [conn] = linkCharges(conn,chemForm,charge)
% Link up species across both ion modes

for n = 1:size(chemForm,1)
    
    fx = strcmp(chemForm,chemForm{n});
    if sum(fx) > 1
        
        % Find the positive ions...
        posI = strcmp(charge(fx),'Positive');
        negI = strcmp(charge(fx),'Negative');
        
        % Determine if we can continue
        if sum(posI) > 0 && sum(negI) > 0
            fx = find(fx);
            posI = fx(posI);
            negI = fx(negI);
            
            for r = 1:numel(posI)                
                for s = 1:numel(negI)
                    
%                     % Check that these values are already 0
%                     vals = [conn(posI(r),negI(s)) conn(negI(s),posI(r))];
%                     if sum(vals) == 100
%                         disp('ok');
%                     elseif sum(vals) > 0
%                         disp('fail');
%                         vals                    
%                     end
                    conn(posI(r),negI(s)) = 30;
                    conn(negI(s),posI(r)) = 30;                    
                end                
            end
            
        end
        %conn(fx,fx) = conn(fx,fx) + 1;
    end
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%