function [factDat,annoInfo] = annotateByGroup(list)
%annotateByGroup - separate variables in this list according to
%differentiating group 

% Indices and unique groups
[unq,~,ind] = unique(list(:,5));
numG = numel(unq);
numV = size(list,1);

% Need to determine the characterisitics for each of the type of lipid
factors = {'Class','Length','Desaturations','Plasmalogen','dCeramide'};
numF = numel(factors);
factDat = cell(numV,numF);

% Loop through each variable
numV = size(list,1);
for n = 1:numV
    
    % Decide on the various parts
    lipid = list{n,6};
    
    % Determine the characteristics. For multiple annotations, some 
    % characteristics may are the same, such as, e.g. PA(O-34:2) and 
    % PA(P-34:1) which are similar in regards to PA and chain length
    [res] = getLipidInfo(lipid);
    
    % So with multiple annotation results, we'll just take the common parts
    % of the annotations, e.g. PA and 34 from above
    for r = 1:size(res,2)
        
        % This is the number/text string
        tmp = res(:,r);
        
        % Compare
        if r == 1 || r == 4 || r == 5
            chk = strcmp(tmp,tmp(1));
        else
            chk = cell2mat(tmp) == cell2mat(tmp(1));
        end
        
        % Then common factor, needs to be saved
        if sum(chk) == numel(tmp)
            
            % Store information that is common to all (or the one)
            % annotation
            factDat(n,r) = {tmp{1}};
        end
    end
end
    
% Want to present one graph for each of the factors, i.e. columns in the
% table. So find all unique values in the full table, then count how many
% in the factDat for each histo group
annoInfo = struct('category',[],'values',[],'histID',[],'table',[]);

% Loop through each of the factors
for n = 1:numF
    
    % Find the non-empty ones
    if n == 2 || n == 3
        fx = cell2mat(factDat(:,n));
        [lev,~,~] = unique(fx);
        fx = ones(size(factDat,1),1);

    else
        fx = ~cellfun(@isempty,factDat(:,n));
        [lev,~,~] = unique(factDat(fx,n));
    end
    numL = numel(lev);
    
    % We want to store the simple count of variables that have each of
    % these characteristics.
    tmp = zeros(numL,numG);
    
    % Loop through each histological group
    for r = 1:numG
        
        % These indices were defined way at the beginning
        fy = ind == r;
        
        % Ensure that we pick non-empty indices from this histo group
        chk = factDat(fx & fy,n);
        
        % For the numeric ones we have to do it a little differently
        if n == 2 || n == 3
            chk = cell2mat(chk);
        end
        
        % 'chk' is a list of the characteristics of these assignments. We
        % find the unique values, and then determine their frequency for
        % the tmp matrix
        [sel,~,idx] = unique(chk);
        numS = numel(sel);
        for s = 1:numS
            
            % Find index of characteristic
            if n == 2 || n == 3
                fz = lev == sel(s);
            else
                fz = strcmp(lev,sel{s});
            end
            
            % Add the number of variables into the matrix
            tmp(fz,r) = sum(idx == s);
        end
        
    end
    
    % Save all of the results together
    annoInfo(n).category = factors{n};
    if n == 2 || n == 3
        annoInfo(n).values = num2cell(lev);
    else
        annoInfo(n).values = lev;
    end
    annoInfo(n).histID = unq;
    annoInfo(n).table = tmp;
    

end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res] = getLipidInfo(lipid)

numA = size(lipid,2);

res = cell(numA,5);

for n = 1:numA
    
    tmp = lipid{1,n};
    
    p1 = strfind(tmp,'(');
    p2 = strfind(tmp,')');
        
    % Class - always the bit before the (
    class = tmp(1:p1-1);
    
    % What is the bit in the bracket?
    chain = ['/' tmp(p1+1:p2-1) '/'];
    
    % Find the slashes to demark the chains
    sl = strfind(chain,'/');
    numSl = numel(sl);
    lnds = zeros(numSl-1,2);
    for r = 1:numSl-1
        
        % Sub-section
        sect = chain(sl(r)+1:sl(r+1)-1);
        
        % Find the colon
        co = strfind(sect,':');
        
        % Length
        lng = sect(1:co-1);
        
        % Check if length entirely numeric
        num = isstrprop(lng,'digit');
        if sum(num) ~= numel(num)        
            
            % Deal with the text bit
            xo1 = lng(~num);
            if strcmpi(class,'cer')
                res{n,5} = xo1;
            elseif strcmp(xo1,'O-')
                res{n,4} = 'O';
            elseif strcmp(xo1,'P-')
                res{n,4} = 'P';
            end
            
            % Deal with the numeric bit
            xo2 = lng(num);
            lnds(r,1) = str2double(xo2);
        else
            lnds(r,1) = str2double(lng);
        end
                
        % Desat
        lnds(r,2) = str2double(sect(co+1:end));        
    end
    
    % Now we use lnds to add to the res cell results
    res{n,1} = class;
    res{n,2} = sum(lnds(:,1));
    res{n,3} = sum(lnds(:,2));
       
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%