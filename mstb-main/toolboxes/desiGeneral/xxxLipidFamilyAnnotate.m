function xxxLipidFamilyAnnotate(~,~,fig,man,mode)
%xxxLipidFamilyAnnotate - annotate spectra and save the families of the
%lipids

% Get the user data
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% What is the ppm tolerance?
ppmTol = str2double(get(man.ppm,'String'));

% This is the database that we want to use
dbPath = deSlash('general/lipid/DBx.mat');

% Use 'mode' to determine the ion mode, which is only relevant when this is
% done for the single ion mode data
if strcmpi(mode,'dual')
    mode = {'pos','neg'};
    %adds = { {'M+H';'M+K'} , {'M-H';'M+Cl'} };
    adds = { {'M+H'} , {'M-H'} };
else
    % Then assume that is 'neg' or 'pos'
    mode = {mode};
    switch mode{1}
        case 'pos'
            adds = {'M+H';'M+K'};
        case 'neg'
            adds = {'M-H';'M+Cl'};
    end
end

% Perform annotation
[db1,ass1] = annotateMZ(dpn.d1.mz,...
    'Polarity',mode{1},...
    'Adducts',adds{1},...
    'ppm',ppmTol,...
    'Database',dbPath);

if isfield(dpn,'d2') % always negative mode
    [db2,ass2] = annotateMZ(dpn.d2.mz,...
        'Polarity',mode{2},...
        'Adducts',adds{2},...
        'ppm',ppmTol,...
        'Database',dbPath);
end

% Now I want to formulate these annotations into something useable in the
% image market...
[fam1] = annoFormat(db1,ass1);
[fam2] = annoFormat(db2,ass2);

% Once we have done enough analyses, then we need to save the results in
% the dpn guidata. 
dpn.d1.anno.db = db1;
dpn.d1.anno.ass = ass1;
dpn.d1.anno.fam = fam1;

if isfield(dpn,'d2')
    dpn.d2.anno.db = []; % don't need to replicate this
    dpn.d2.anno.ass = ass2;
    dpn.d2.anno.fam = fam2;
end

% Update the guidata
guidata(fig.fig,dpn);

% Now we could update the popup menu in the toolbox
if isfield(dpn,'d2')
    allFam = unique([fam1(:); fam2(:)]);
else
    allFam = unique(fam1(:));
end
set(man.family,'String',allFam,'Value',1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fam] = annoFormat(db,ass)
% Make the annotations useful in terms of lipid families

fam = repmat({'None'},[size(ass.match,1),3]);

% Loop through each one
for n = 1:size(fam,1)
    
    % Skip empty ones
    if isempty(ass.match{n})
        continue;
    end
    
    % Indices in the DB
    idx = ass.match{n};
        
    % Determine the class(es) that this feature could be, and only include
    % if all classifications are the same!
    c1 = unique(db.Class1(idx));
    if numel(c1) == 1
        fam(n,1) = c1;
    end
    
    c2 = unique(db.Class2(idx));
    if numel(c2) == 1
        fam(n,2) = c2;
    end
    
    c3 = unique(db.Class3(idx));
    if numel(c3) == 1
        fam(n,3) = c3;
    end
    
end
    



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
