function [ meta ] = mtspMVA(res)
%mtspMVA - prepare the variables based on annotations

% Determine the unique annotations...
[all] = filterAnnotations(res,'neg');

% Make the matrix of variables
[data,meta] = makeMatrix(all);

figure; imagesc(data);

doPCA(data,meta)


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [all] = filterAnnotations(res,charge)

% Combine res into a single large table
all = vertcat(res{:});

% Filter out MRMs which are false in the FDR
filt = strcmp(all(:,10),'True');
all = all(filt,:);

% What about MSM > 0.9
msm = cell2mat(all(:,8));
filt = msm > 0.9;
all = all(filt,:);

% Include only annotations from the same DB - either HMDB or ChEBI
filt = strcmp(all(:,1),'HMDB');
all = all(filt,:);

% Decide which charged species to keep
switch charge
    
    case 'pos'
        fx = strcmp(all(:,12),'Positive');
        all = all(fx,:);
    case 'neg'
        fx = strcmp(all(:,12),'Negative');
        all = all(fx,:);
        
    otherwise
        % Do nothing - use both
end
        
        
    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [matr,meta] = makeMatrix(all)

% Now determine the unique annotations - let's use the adduct as well
[unqA,firstA,~] = unique(all(:,4));
numA = numel(unqA);

% Now determine the unique files
[unqF,~,~] = unique(all(:,2));
numF = numel(unqF);

% Create an empty matrix
matr = zeros(numF,numA);

% Loop throuhgh and populate
for n = 1:size(all,1)
    
    iF = strcmp(unqF,all(n,2));
    iA = strcmp(unqA,all(n,4));
    
    matr(iF,iA) = 1;
    
end

% Tissue types
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
meta.tissue = tissue;

% Variable m/z value
meta.mz = cell2mat(all(firstA,5));

meta.anno = unqA;
meta.file = unqF;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doPCA(data,meta)

[ll,ss,ee] = pca(data);
ee = 100 * ee / sum(ee);


% Plot observations
scatterPlot(ss(:,1:2),meta.tissue,'','PC1','PC2');

figure; scatter(ll(:,1),ll(:,2),80,meta.mz,'o','filled');
cb = colorbar;

figure; stem(meta.mz,ll(:,1));


% MMC classification
[a,b,c,d,e] = kfoldMMC(data,meta.tissue,10);
makeConfMat([],b,parula(5),[],d)
scatterPlot(a(:,1:2),meta.tissue,'','LV1','LV2');



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
