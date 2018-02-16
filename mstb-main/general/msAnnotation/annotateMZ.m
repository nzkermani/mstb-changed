function [lm,ass] = annotateMZ(mz,varargin)
% annotateMZ - an improvement on peakAnnotationBatch

tic

% Get the varargin
[opts] = readArgsData(varargin);

% Read in the database
if ~exist(opts.database,'file')
    disp('Cannot locate database');
    disp('Are you in the correct folder');
    return
end
tmp = open(opts.database);
lm = tmp.lm;
clear tmp

% Depending on the adducts selected by the user, we need to calculate the
% relevant m/z values from the masses
[add] = determineAdducts(opts);

% Convert the masses to relevant m/z values
[mzgrid] = mass2mz(lm.Mass,add);
[~,numA] = size(mzgrid);

% Calculate the low and high m/z values of the list based on ppm tolerance
tols = opts.ppmTol .* mz ./ 1e6;

% For each ion in the list provided, we need to match it according to the
% m/z values derived from the database
numV = numel(mz);

% Vector to save the matches
match = cell(numV,numA);

% What about saving information regarding which entries in the database are
% annotated. Luisa, for example, would like to see how lipid classes change
% over time
numDB = size(lm.Formula,1);
dbMatch = false(numDB,1);

% Loop through each variable
for n = 1:numV
    
    % See if any DB entries match mz(n)
    tmp = abs(mzgrid - mz(n));
    tmp = tmp < tols(n);
    
    % Find matches and save their indices
    [fx,fy] = find(tmp == 1);
    if numel(fx) == 0
        continue;
    end

    % Loop through poss
    for r = 1:numA
        fz = fy == r;
        
        if sum(fz) > 0
            match{n,r} = fx(fz);
        end
    end
    
    % Save to the dbMatch list
    dbMatch(fx,1) = true;
end
   
% What kind of output do we want?
ass.match = match;
ass.adduct = add;
ass.opts = opts;
ass.dbMatch = dbMatch;

toc

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = readArgsData(argsin)
% Read the arguments and then the data if it wasn't passed

% Define the defaults here
opts.polarity   = 'neg';
opts.adduct     = {'M-H','M+Cl'};
opts.ppmTol     = 5;
opts.database   = deSlash('general/lipid/DB1.mat');

% Run through each pair
nArgs = length(argsin);
for i = 1:2:nArgs
    if strcmpi('adduct',argsin{i}) || strcmpi('adducts',argsin{i})
        tmp = argsin{i+1};
        if iscell(tmp)
            opts.adduct = tmp;
        end
        
    elseif strcmpi('polarity',argsin{i})
        tmp = argsin{i+1};
        if strcmp(tmp(1),'p')
            opts.polarity = 'pos';
        elseif strcmpi(tmp(1),'n')
            opts.polarity = 'neg';
        end
        
    elseif strcmpi('tolerance',argsin{i}) || strcmpi('ppm',argsin{i})
        tmp = argsin{i+1};
        if isnumeric(tmp)
            opts.ppmTol = tmp;
        end
        
    elseif strcmpi('db',argsin{i}) || strcmpi('database',argsin{i})
        tmp = argsin{i+1};
        
        if strcmpi(tmp,'dipa') || strcmpi(tmp,deSlash('general/lipid/DB1.mat'))
            opts.database = deSlash('general/lipid/DB1.mat');
        elseif strcmpi(tmp,'luisa') || strcmpi(tmp,deSlash('general/lipid/DBx.mat'))
            opts.database = deSlash('general/lipid/DBx.mat');
        else
            error('Unknown database issue');
        end
        
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [add] = determineAdducts(opts)

% Get the full adduct list
[list,name,elem] = adductLists(opts.polarity(1));

% Now we need to match name with opts.adduct
numA = numel(opts.adduct);
idx = zeros(numA,1);

for n = 1:numA    
    fx = strcmp(name,opts.adduct{n});    
    idx(n,1) = find(fx == 1);    
end

% Trim the list
add.list = list(idx,:);
add.name = name(idx,:);
add.elem = elem(idx,:);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function junk

% Calculate m/z tolerances for each database m/z entry
numV = numel(mz);

% Empty places for storing the results: keep: mzVar / dbMatch
initMatch   = 2000;
cmpList     = zeros(initMatch,2);
varAnno     = zeros(numV,1);
annoNam     = cell(numV,1);
annoCls     = cell(numV,1);

% These are places for the storage of isotopic information
isIso       = zeros(numV,2); % if true contains m/z and index of M
iso1        = zeros(numV,2); % mz and index of M+1
iso2        = zeros(numV,2); % mz and index of M+2

% We have a series of mz values from the experiment (mz) and a list of
% database entries with mz values (dbMZ) and names (dbNam). The task is to
% match the two to each other...
tic
numMatch = 0;
for n = 1:numV
    %disp([int2str(n) '---' num2str(mz(n))]);
    
    % See if any DB entries match mz(n)
    tmp = abs(db.mz - mz(n));
    tmp = tmp < tols;
    
    % If there are any hits
    if sum(tmp) > 0
        
        % Get their indices
        [~,fy] = find(tmp > 0);
        
        for r = 1:numel(fy)
            numMatch = numMatch + 1;
            cmpList(numMatch,:) = [n fy(r)];% mz(n) dbMZ(fy(r))];
            
            if r == 1
                annoNam{n,1} = char(db.nam(fy(r)));
                annoCls{n,1} = [char(db.c1(fy(r))) ',' char(db.c2(fy(r)))];
            else
                annoNam{n,1} = [char(annoNam{n,1}) ' / ' char(db.nam(fy(r)))];
                annoCls{n,1} = [char(annoCls{n,1}) ' / ' char(db.c1(fy(r))) ',' char(db.c2(fy(r)))];
            end
        end
        
        % This is indended as a vector which can be used to identify how 
        % many annotations each variable has.
        varAnno(n,1) = varAnno(n,1) + numel(fy);
    else
        annoNam{n,1} = '';
    end
    
    % Once it has done the main matching algorithm, then we need to get the
    % matching of isotopic peaks done.  Of course a peak can have both an
    % assignment in the DB and be considered as an M+n peak, but somewhere
    % (later) their dual statuses will have to be resolved.
    
    % Only proceed if this peak hasn't been nabbed as an M+[0.5/1/2] peak
    if isIso(n,2) == 0
    
        tol = mztolDB * mz(n) / 1e6;

        % These are the margins from the mz(n) peak being investigated
        isol = (mz - tol) - mz(n);
        isoh = (mz + tol) - mz(n);

        isov = [isol' isoh'];

        % Find local matches
        [fHalf,~] = find(isov(:,2) > cdiff/2 & isov(:,1) < cdiff/2);
        [fOne, ~] = find(isov(:,2) > cdiff*1 & isov(:,1) < cdiff*1);
        [fTwo ,~] = find(isov(:,2) > cdiff*2 & isov(:,1) < cdiff*2);

        % Start with halves as these will be the first, then move onto 1,
        % but leave 2 if there is no one
        if ~isempty(fHalf)
                        
            % Look for + 0.5
            if numel(fHalf) == 1
                val = 1;                
            else
                mn = mean(isov(fHalf,:),2) - cdiff/2;
                [~,val] = min(abs(mn));                
            end     
            
            if ~isempty(fHalf)
                % Add into the matrix
                iso1(n,:) = [mz(fHalf(val)) fHalf(val)];

                % Now we need to 'tell' the isIso matrix that this is an
                % isotope peak so that we ignore it when we come to it...
                isIso(fHalf(val),:) = [mz(n) n];
            end            
            
            % Now look for +1
            if numel(fOne) == 1
                val = 1;                
            else
                mn = mean(isov(fOne,:),2) - cdiff;
                [~,val] = min(abs(mn));                
            end      
            
            if ~isempty(fOne)
                % Add into the matrix
                iso2(n,:) = [mz(fOne(val)) fOne(val)];

                % Now we need to 'tell' the isIso matrix that this is an
                % isotope peak so that we ignore it when we come to it...
                isIso(fOne(val),:) = [mz(n) n];
            end
            

        elseif ~isempty(fOne)

            % Then can focus on the +1 +2 combination
            if numel(fOne) == 1 
                val = 1;                
            else
                mn = mean(isov(fOne,:),2) - cdiff;
                [~,val] = min(abs(mn));                
            end
            
            if ~isempty(fOne)
                % Add into the matrix
                iso1(n,:) = [mz(fOne(val)) fOne(val)];

                % Now we need to 'tell' the isIso matrix that this is an
                % isotope peak so that we ignore it when we come to it...
                isIso(fOne(val),:) = [mz(n) n];
            end
            
            
            % Now need to look for M+2 peaks - same procedure as above
            if numel(fTwo) == 1
                val = 1;
            else
                mn = mean(isov(fTwo,:),2) - 2 * cdiff;            
                [~,val] = min(abs(mn));                
            end
            
            if ~isempty(fTwo)
                % Add into the matrix
                iso2(n,:) = [mz(fTwo(val)) fTwo(val)];
            
                % Now we need to 'tell' the isIso matrix that this is an
                % isotope peak so that we ignore it when we come to it...
                isIso(fTwo(val),:) = [mz(n) n];
            end

        elseif ~isempty(fTwo)
            % Do nothing here, as expect/demand to see an M+1 peak to even
            % consider seeing the M+2 peak. This if this is the only match then
            % it is as good as ignored.
        end
    end    
    
    
    
    
    
    
end
disp(['Number of matches: ' int2str(numMatch) char(9) 't = ' num2str(toc) 's']);

% Trim to remove unused rows
cmpList = cmpList(1:numMatch,:);

% So now we have a short list of the matches between the actual spectral
% variables (col1) and the db entires (col2). Any other information can be
% scrubbed from the appropriate place, e.g. the intensity matrix.

% Format the output into a structure and then return...
ass.list    = cmpList;  % matches between mzVar (:,1) and db entry (:,2)
ass.db      = db;       % db structure from dbRead
ass.ppmTol  = mztolDB;  % provided mz tol
ass.annoNum = varAnno;  % how many matches to each variable
ass.annoNam = annoNam;  % string containing match(es) of annotations
ass.annoMZ  = mz;       % the mz values provided (not sure necessary?)
ass.annoCls = annoCls;  % the class/subclasses of the annotations

% Isotopic information in here
ass.isIso = isIso;
ass.iso1  = iso1;
ass.iso2  = iso2;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
