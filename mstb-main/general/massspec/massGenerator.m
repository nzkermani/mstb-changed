function [ mz ] = massGenerator( varargin )
%massGenerator - for a series of stuff, generate the monoisotopic mass
%values for it. This will be a very good function for generating the
%possible fatty acids, e.g. R+CxHy with increasing length and
%desaturations.
%
% The following elements are supported. This is because of the isotopic
% calculations which follow it. If you change the element, then there can
% be no adequate isotopic stuff that follows... (Unless I make it all a
% single function, thus cementing the change to a single place! To Do.
% C H Li N O F Na Mg Si P S Cl K Ca As Se Br I

% James S. McKenzie, Imperial College, London. 2014.

% First get the inputs sorted out...
[opts] = getVarArgin(varargin);

% There are various ways to use this function. Simplest of all is to
% provide a single formula of something, and we'll calculate the mass of it
% using accurate monoisotopic masses
if ~isempty(opts.formula)
    elem = parseFormula(opts.formula);    
    elem   = mzCalculate(elem);
    [isoDist] = isotopicCalculator(elem.qty)
    
    % Convert a single accurate monomass to a series of length(isoDist) by
    % simply adding one each time...
    isoMZ = ones(1,numel(isoDist)) * elem.mz + [0:1.0034:numel(isoDist)-1];
    
    figure; stem(isoMZ,isoDist);
    mz = elem.mz;
    return
end

% How about a geometric like progression of alkanes from CH3(CH2)n- for a
% range of n. Then can desaturate by removing H2
if ~isempty(opts.carbons)
    [allE] = alkaneSeries(opts);
    
    % And if there are some desaturations to be included?
    if ~isempty(opts.desats)
        
        % Then progress the alkane series by performing as many
        % desaturations as possible...
        [mz] = alkaneDesaturase(opts,allE);
        
        % Output the file...
        tt = cell2table(mz(:,1:4));
        writetable(tt,'tmp.txt','delimiter','\t',...
            'WriteVariableNames',0)
        
        % Make an excellent plot of lipid isotopic distrubtions
        distPlot(mz,opts);
        return
    else
        return
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = getVarArgin(argsin)

opts.formula    = [];
opts.head       = 'H';
opts.carbons    = [];
opts.desats     = [];
opts.ions       = 'M';
opts.style      = 'M';

% Now run through the list of arguments
nArgs = length(argsin);
for i=1:2:nArgs
    if strcmpi('Head',argsin{i})
        opts.head   = argsin{i+1};
        
    elseif strcmpi('Carbons',argsin{i})
        opts.carbons = argsin{i+1};
        if numel(opts.carbons) == 1
            opts.carbons = [opts.carbons opts.carbons];
        end
        
    elseif strcmpi('Desats',argsin{i})
        opts.desats = max(argsin{i+1});

    elseif strcmpi('Formula',argsin{i})
        opts.formula = argsin{i+1};

    elseif strcmpi('Ion',argsin{i})
        opts.ions = argsin{i+1};
        
    elseif strcmpi('Style',argsin{i})
        opts.style = argsin{i+1};

    
    end        
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [elem] = defineElements
% Gather the elementarily necessary information for the function. In simple
% form this is the list of relevant elements: CHNOPS

elem.symb = {'C','H','Li','N','O',...
    'F','Na','Mg','Si','P',...
    'S','Cl','K','Ca','As',...
    'Se','Br','I'};
elem.mono = [12,...     % C12
    1.00782503027,...   % H1
    6.015122795,...     % Li3 - this is only 8% abundant
    14.0030740048,...   % N14
    15.99491461956,...  % O16
    18.99840322,...     % F19
    22.9897692806,...   % Na23
    23.985041700,...    % Mg25
    27.9769265325,...   % Si28
    30.9737663,...      % P31
    31.97207100,...     % S32
    34.96885268,...     % Cl35
    38.96370668,...     % K39
    39.96259098,...     % Ca40
    74.9215965,...      % As75
    73.9224764,...      % Se74 - 0.8% abundant
    78.9183371,...      % Br79
    126.904473];        % I127

numE = numel(elem.mono);
elem.qty = zeros(1,numE);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [elem] = parseFormula(form)
% Convert a string of chemical elements into a single homogeneous string
% that you can use in the function. Will require a defined list of
% elements, in order...

% Define the elements
[elem] = defineElements;

% Shove a 1 on the end if absent...
if isstrprop(form(end),'lower') || isstrprop(form(end),'upper')
    form(end+1) = '1';
end

% Calculate the population: text / numeric
txt = isstrprop(form,'upper');
low = isstrprop(form,'lower');
num = isstrprop(form,'digit');


% How many chars in the formula?
lng = length(form);

% How many element entries are there?
numE = sum(txt);
ex = find(txt == 1);

for n = 1:numE
    
    % This is the element
    st = ex(n);    
    if low(ex(n)+1) == 1
        fn = ex(n)+1;
    else
        fn = st;
    end
    ee = form(st:fn);
    
    % Now the numbers... What is beyond fn?
    if num(fn+1) == 0
        qq = 1;
    else
        jj = fn+1;
        while num(jj) == 1 && jj < lng %|| jj <= lng-1
            jj = jj + 1;
        end
        
        if jj ~= lng
            jj = jj - 1;
        end
        qq = str2double(form(fn+1:jj));
    end
    
    % Now need to put into the qty vector in elem.qty
    ei = strcmp(elem.symb,ee);
    if sum(ei) == 0
        error('Invalid element - don''t know what to do with it');
    end
    elem.qty(ei) = elem.qty(ei) + qq;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [elem] = mzCalculate(elem)
% Determine the m/z of the parsed formula

indiv   = elem.mono .* elem.qty;
elem.mz = sum(indiv);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [allE] = alkaneSeries(opts)
% A geometric progression of increasing chain length in alkanes...

% The 'formula' for linear alkanes of form R-H is CH3(CH2)x-1, which
% translates to be C_x H_(2x+1)
n = 0;
allE = struct('C',[],'elem',[]);
for i = min(opts.carbons):max(opts.carbons)
    n = n + 1;
    
    allE(n).elem = defineElements;
    allE(n).C    = i;
    
    % Add in carbons...and hydrogens
    allE(n).elem.qty(1) = i;
    allE(n).elem.qty(2) = i*2 + 1;
    
    % Calculate the mass
    allE(n).elem = mzCalculate(allE(n).elem);
        
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [op] = alkaneDesaturase(opts,allE)
% Like the enzymes, let's take an alkane series and sequentially desaturate
% them to produce all possible alkenes for each chain length...

% Also want to make a cell structure like thing which can be exported as a
% table in Excel and used like the current database.

% How many desats?
numD = max(opts.desats);

listD = 0:1:max(opts.desats);
listL = min(opts.carbons):max(opts.carbons);

rn = cell(numel(listL),1);
vn = cell(numel(listD),1);

% Grid...
grid = zeros(size(allE,2),numD+1);

% Table thing for export of all possibilities...
tabs = zeros(numel(grid),numel(allE(1).elem.qty)+1);


% What is the mass of 2H?
tmp = defineElements;
dst = 2 * tmp.mono(2);
k = 0;

% THis is the head group which we should include...
head  = parseFormula(opts.head);
head  = mzCalculate(head);

subt = zeros(size(head.qty));
subt(2) = 2;

for n = 1:size(allE,2)
    
    % Max permissible desaturations?
    if listL(n) == 1
        maxD = 0;
    else
        maxD = min([opts.desats floor(allE(n).elem.qty(1) / 2)]);
    end
    
    % Add the saturated form in first
    k = k + 1;
    grid(n,1) = allE(n).elem.mz;
    rn{n,1} = ['C_' int2str(listL(n))];
    tabs(k,1:end-1) = allE(n).elem.qty;
    
    % Now for sequential desaturations
    for r = 1:maxD
        if r == 1
            grid(n,r+1) = allE(n).elem.mz - dst;
            
        else
            grid(n,r+1) = grid(n,r) - dst;
        end
        
        % All this should be entered into the table as well...
        k = k + 1;
        tabs(k,1:end-1) = allE(n).elem.qty - (r*subt);
        tabs(k,end)     = r;
        
    end
    
        
    
end

% These are the column labels showing desaturations
for r = 1:numel(listD)
    vn{r,1} = ['d' int2str(listD(r))];
end

% Trim the table to remove additional entries
tabs = tabs(1:k,:);

% Here is where we incorporate the head group into the figures, as
% previously it had been omitted to make for an easier alkane series
inds  = grid ~= 0;
grid2 = grid + head.mz;
grid2 = grid2 .* inds;
tabs  = bsxfun(@plus,tabs,[head.qty 0]);

frags = array2table(grid,'VariableNames',vn,'RowNames',rn)
full  = array2table(grid2,'VariableNames',vn,'RowNames',rn)

% What is the ion form? We'll count the electron's mass as neglibible
ionForm = defineElements;
switch opts.ions
    case {'[M-H]'; '[M-H]-'}
        ionForm.qty(2) = -1;        
        
    case {'[M+FA-H]'; '[M+FA-H]-'}
        ionForm.qty(1) = 1;
        ionForm.qty(2) = 1;
        ionForm.qty(5) = 2;
        
    otherwise
        error('Incorrect ion specified');
end
        
        

% Now with the table of elements in each of the permutations, let's make a
% good table with masses, names, and actual IONs, using the Ion input...
numI = size(tabs,1);
genr = defineElements;
op = cell(numI,6);
for n = 1:numI
    
    % Take the entry from the table and add to the object thing
    genr.qty = tabs(n,1:end-1);
    
    % Now put in the ion / loss or gain...
    genr.qty = genr.qty + ionForm.qty;  
    
    % Here calculate the mz of the ion
    genr = mzCalculate(genr);
    
    % Generate a text chem form
    frm = elem2form(genr);
    
    % Now calculate the isotopic distribution of the species...
    isodist = isotopicCalculator(genr.qty);
    isomzs  = ones(1,numel(isodist)) * genr.mz + [0:1.0034:numel(isodist)];
    
    op{n,1} = genr.mz;    
    op{n,2} = [opts.style '(' int2str(genr.qty(1)-ionForm.qty(1)) ':' int2str(tabs(n,end)) ')'];
    op{n,3} = frm;
    op{n,4} = opts.ions;
    
    op{n,5} = isomzs;
    op{n,6} = isodist;
    
    
end
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [form] = elem2form(elem)
% Convert to a chemical formula

numE = numel(elem.qty);
for n = 1:numE    
    if elem.qty(n) ~= 0
        if elem.qty(n) == 1
            fm = elem.symb{n};
        else
            fm = [elem.symb{n} int2str(elem.qty(n))];
        end
        
        if ~exist('form','var')
            form = fm;
        else
            form = [form fm];
        end
    end
end
        

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function distPlot(mz,opts)
% Make a plot of the isotopic distributions...

cols = jet(opts.desats+1)

numL = size(mz,1)

figure; hold on;

for n = 1:numL
    
    % This is the colour index
    ci = mod(n-1,size(cols,1))+1;
    
    % How many carbons / desats
    cs = strfind(mz{n,2},'(')
    cc = strfind(mz{n,2},':')
    cf = strfind(mz{n,2},')')
    
    nc = str2double(mz{n,2}(cs(1)+1:cc(1)-1))
    nd = str2double(mz{n,2}(cc(1)+1:cf(1)-1))
    
    if nd == 0
        marker = 'o';
        sz = 15;
    else
        marker = 's';
        sz = 10;
    end
    
    
    
    
    stem(mz{n,5},mz{n,6},'Marker',marker,...
        'MarkerFaceColor',cols(nd+1,:),...
        'MarkerEdgeColor',cols(nd+1,:),...
        'MarkerSize',sz);
    
    
    
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%