function [ass] = peakAnnotationBatch(mz,mztolDB,dbP)
% peakAnnotationBatch - an improvement on Ottmar's version, as this is
% considerably quicker. The output format has changed somewhat. The
% principle is based on the simple comparison of m/z values that agree
% within a ±ppm window. There may be multiple matches. The various elements
% of the output structure ass(igned) are explained below.
%
%   mz:         matrix mz values
%---X:          matrix of peak intensities of mz
%   mztolDB:    mz tolerance for peak annotation in ppm
%   DBPath:     path to database excel file containing mz values in first
%               column
%
%   ass(ign)    structure, for detailed info of each element see the
%               desciptions provided below

% Simple defaults here
if nargin == 2
    
    %if max(mz) <= 400
    %    dbP = ['sourcePackages' filesep 'FADatabaseV02.xlsx'];        
    %else
        dbP = ['packages' filesep 'msDatabases' filesep 'LipidDatabaseV03.xlsx'];
        
        %warning('Changed the default DATABASE');
        %dbP = ['sourcePackages' filesep 'NimaFattyHumans.xlsx'];
        
    %end
end

% Need to add some work on the isotopic procedure. Will look at peaks
% doubly or singly charged, and will actively look for the M+1 and M+2
% peaks. Predominantly expecting z = 1, but doubly charged z = 2 will be
% accommodated. The likelihood of finding an isotopic peak depends on the
% accuracy of the instrument.
cdiff = 1.00335;

% Read in the database
[db] = dbRead(dbP);

% Calculate m/z tolerances for each database m/z entry
tols = mztolDB .* db.mz ./ 1e6;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [db] = dbRead(dbP)
% Read the excel spreadsheet masquerading as the database. Get the useful
% bits from it, and save to a structure.

% Read in from excel
[~,~,edb] = xlsread(dbP);

db.mz    = cell2mat(edb(:,1)'); % m/z values of ions
db.nam   = edb(:,2)';           % name of entries
db.form  = edb(:,3)';           % their formulae
db.ion   = edb(:,4)';           % and their ion forms

% Subsequent columns are classes and subclasses and subsubclasses...
numC = size(edb,2);
if numC > 4
    for n = 5:6
        tmpName = ['c' int2str(n-4)];
        db.(tmpName) = edb(:,n)';
    end
end
db.numC = numC-4;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

