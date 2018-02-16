function [gr] = mtspGraph(res)
%mtspGraph - make a graph of some type to show connectivity of samples by
%molecular annotations

% Combine res into a single large table and then determine connectivities
all = vertcat(res{:});

% Filter out MRMs which are false in the FDR
filt = strcmp(all(:,10),'True');
all = all(filt,:);

% What about MSM > 0.9
msm = cell2mat(all(:,8));
filt = msm > 0.95;
all = all(filt,:);

% Include only annotations from the same DB - either HMDB or ChEBI
filt = strcmp(all(:,1),'HMDB');
all = all(filt,:);

% We need to add a column for the charge of the ion - in most cases this is
% negative, although there are a few positive samples (ovarian)


% Annotation graph
%[gr] = fileGraph(all);

[gr] = annoGraph(all);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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