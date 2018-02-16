function [ all ] = hnak(op,ass,all)
%hnak - m/z analysis for the peaks found with H/Na/K adducts. For Luisa.
%I'll put a better description in here when I've worked out how the
%function is likely to work.
%
% INPUTs (* denotes optional input)
% op    - from lipid analysis toolbox
% ass   - also from lipid analysis toolbox
% all*  - the results structure, for launching the GUI again
%
% OUTPUTs
% all   - structure containing all of the results, for reloaded the GUI
%
% FILEs
% txt   - a txt file appended with the date/time is output in the current 
%         folder
%
% CHANGE LOG
%
% 06/05/15
% Lipid annotations are now exported; these are matches in the database
% to the M+H adduct.
% Boxplots are now drawn in a GUI for each lipid / species.
% Can save the results, and reload for visualisation.
%
% 07/05/15
% Now boxplots differentiate between histological groups too.
%
% James McKenzie, 2015.
%
%
% TO DO
% - Incorporate the 'minFind' value, so as to not exclude variables where
%   only some adducts are found.


% The most important part is the specification of which adducts are to be
% searched for.  Note that the differences are relative to the first 
% specified ion rather than to the actual mass.
add.name    = {'H','Na','K'};
add.diff    = [0 21.981944 37.955882];
add.thresh  = 5;
add.minFind = 2;


% Find interesting ions, with indices of H/Na/K ions
if nargin == 2
    [res] = mzSearch(op,add);

    % Now we need to run through the samples and find the ratios of these
    [tab.dat,tab.hds] = dataScrub(op,add,res,ass);

    % Let's output the full table and see what Luisa wants to do with it...
    tabOP(tab);

    % How about a structure to save everything...
    all.res = res;
    all.tab = tab;
    all.add = add;
   
end    

% How about a gui-like window for the visualisation of the results, such as
% drawing some box plots...
hnakGUI(all);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res] = mzSearch(op,add)
% Run through the cmz vector and identify peaks with the desired pattern

% How many adducts?
numA = numel(add.diff);

% How many variables?
numV = numel(op.cmz);

% Somewhere to store the results - indices of matching peaks!
res = NaN(numV,numA);

% Loop through each variable
for n = 1:numV
    
    % Replicate the mz vector
    tmpMZ = repmat(op.cmz',1,numA);

    % Subtract the current m/z variable
    tmpMZ = tmpMZ - op.cmz(n);
    
    % Now remove the adduct differences
    tmpMZ = bsxfun(@minus,tmpMZ,add.diff);
    
    % Convert these into ppm values for a more instructive tolerance
    % comparison
    ppm = abs(1e6 * bsxfun(@rdivide,tmpMZ,op.cmz'));
    
    % Set larger ppm vales to NaN to make it easy to spot if there are no
    % potential matches
    ppm(ppm > add.thresh) = NaN;
    
    % Find the minimum values and indices for peaks...
    [ppm,ind] = min(ppm,[],1);
    
    % Now save the results in a matrix
    ind(isnan(ppm)) = NaN;
    res(n,:) = ind;
    
end

% Let's filter out the non-results, i.e. we want only peaks which have all
% of the adducts specified in the list
nnn = sum(isnan(res),2);

% Which ones have no nans in?
fx = nnn == 0;

% Trim the list!
res = res(fx,:);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tab,hds2] = dataScrub(op,add,res,ass)
% Run through the actual data, getting the suitable data and make it
% presentable

numObs = size(op.histID,1);
[numIon,numAdd] = size(res);

[unqH,~,~] = unique(op.histID);
numHist = numel(unqH);


% Prepare the column headings...
init = {'Patient ID','Sample ID','Hist ID',['Anno: M*' add.name{1}]};
t = 0;
for r = 1:numel(add.name)
    for s = r+1:numel(add.name)
        t = t + 1;
        gener{1,t} = [add.name{r} '/' add.name{s}];
    end    
    gener{1,t+1} = add.name{1};
end

% Concatenate the headings
hds2 = {init{:} add.name{:} add.name{:} gener{:}};

% Now that we have hacked the headings together, we can start to form the
% empty table, and then fill it...
tab = cell(numObs*numIon,numel(hds2));
n = 0;

% Loop through all of the observations...
for o = 1:numObs
    
    % Add then for each group of ions
    for i = 1:numIon
        r = 1;
        n = n + 1;
        
        % patient ID
        tab{n,r} = op.pID{o};
        r = r + 1;
        
        % sample ID
        tab{n,r} = op.sampleID{o};
        r = r + 1;
        
        % hist ID
        tab{n,r} = op.histID{o};
        r = r + 1;
        
        % Now how about squeezing in any annotations that may exist - there
        % may be more than one
        annoInd = res(i,1); % this is the suspected H adduct        
        annoTxt = ass.annoNam{annoInd};
        tab{n,r} = annoTxt;
        r = r + 1;
        
        % Ion indices
        for q = 1:numAdd
            tab{n,r+q-1} = res(i,q);
        end
        r = r + numAdd;
        
        % Ion m/z values
        for q = 1:numAdd
            tab{n,r+q-1} = op.cmz(res(i,q));
        end
        r = r + numAdd;
        
        % Now the values (calculate them first)
        rats = cell(1,numel(gener));
        u = 0;
        for s = 1:numAdd
            for t = s+1:numAdd
                u = u + 1;             
                rats{1,u} = op.XPeaks(o,res(i,s)) / op.XPeaks(o,res(i,t));
            end
        end
        rats{1,u+1} = op.XPeaks(o,res(i,1));
        
        % Add them in...
        for q = 1:numel(gener)
            tab{n,r+q-1} = rats{1,q};
        end
        %r = r + numel(gener);        
    end
end
       



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tabOP(tab)
% Output the table in a text format, which can be opened with Excel

fNam = [pwd filesep 'Luisa-Lipid-' datestr(now,'yymmdd-HHMMSS') '.txt'];

fID = fopen(fNam,'w');

format = {'%d\t','%0.4f\t'};

% Print the column headings
[numRow,numCol] = size(tab.dat);
for n = 1:numCol
    fprintf(fID,'%s\t',tab.hds{n});
end
fprintf(fID,'\n');

% Now the main bulk of the table...
for r = 1:numRow
    
    fprintf(fID,'%s\t%s\t%s\t%s\t',...
        ['[' tab.dat{r,1} ']'],tab.dat{r,2},tab.dat{r,3},tab.dat{r,4});

    for c = 5:numCol        
        
        if c <= 6
            tf = format{1};
        else
            tf = format{2};
        end
        
        if isnan(tab.dat{r,c}) || isinf(tab.dat{r,c})
            fprintf(fID,'\t');
        else
            fprintf(fID,tf,tab.dat{r,c});
        end
    end
    fprintf(fID,'\n');
end

fclose(fID);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hnakGUI(all)
% Draw a gui for showing the results of the HNAK analysis. All we need at
% the beginning is a table, and a drop-down for the showing of lipid box
% plots or something

% Determine the unique lipid annotations
isEmp   = ~cellfun(@isempty,all.tab.dat(:,4));
actLip  = all.tab.dat(isEmp,4);
unqLip  = unique(actLip);

% Also want the indices of the unique species (annotation or not)
allInd = cell2mat(all.tab.dat(:,5));
unqInd = unique(allInd);

% Draw a window
[fig] = drawGUI(unqLip,unqInd);

% Add the information into the table...
set(fig.tab,'Data',all.tab.dat,...
    'ColumnName',all.tab.hds);

% Add the callback functionality
set(fig.dd1,'Callback',{@ddcb,fig,all});
set(fig.dd2,'Callback',{@ddcb,fig,all});

% Initialisation
ddcb(fig.dd1,[],fig,all);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = drawGUI(list,inds)
% Draw the GUI...

f0 = findobj('Tag','hnak');
close(f0);

fig.fig = figure('Name','HNaK Analysis',...
    'Tag','hnak',...
    'Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5],...
    'Menubar','none',...
    'Number','off',...
    'Toolbar','figure');

fig.tab = uitable('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.05 0.1 0.4 0.8]);

fig.ax = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.55 0.1 0.4 0.8]);

fig.dd1 = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.55 0.925 0.2 0.05],...
    'Style','popupmenu',...
    'String',cat(1,'All',list),...
    'Value',1,...
    'Tag','dd');

fig.dd2 = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.75 0.925 0.2 0.05],...
    'Style','popupmenu',...
    'String',int2str([0; inds]),...
    'Value',1,...
    'Tag','dd');


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ddcb(src,~,fig,all)
% Dropdown callback function - this is what happens when you click the
% button...

% Find the other dd's handle, and set it's value back to 1...
f0 = findobj('Tag','dd');
delF = setdiff(f0,src);
set(delF,'Value',1);


% Get the selection...
lab = get(src,'String');
ind = get(src,'Value');

cla(fig.ax,'reset');

if ind == 1
    % This is 'all'
    
    set(fig.ax,'Visible','off');
    set(fig.tab,'Data',all.tab.dat);
    return
end

% Otherwise, we can continue... Need to modify to allow for either box
if iscell(lab)
    
    lipid = lab(ind);
    
    % Find this specific lipid in the table of data
    inds = strcmp(all.tab.dat(:,4),lipid);

elseif ischar(lab)
    lipid = str2num(lab(ind,:));
    
    inds = cell2mat(all.tab.dat(:,5)) == lipid;
    
end


% Trim the dataset
data = all.tab.dat(inds,:);

% Update in the table
set(fig.tab,'Data',data);

% Now we need to get some stuff out of the data table...
sfind = strfind(all.tab.hds,'H/');
sfind = ~cellfun(@isempty,sfind);
sfind = find(sfind == 1);
inds = [sfind(1) numel(all.tab.hds)-1];

% How big is the data?
numO = size(data,1);
numA = inds(2)-inds(1)+1;
bpd = NaN(numO*numA,1);
bpi = cell(size(bpd,1),2);

% Now loop through
s = 1;
for n = inds(1):inds(2)
    
    % High index
    t = s + numO - 1;
    
    % Add label to bpi
    bpi(s:t,1) = all.tab.hds(n);
    
    % And histID
    bpi(s:t,2) = data(:,3);
    
    % Now plonk in the data    
    bpd(s:t,1) = cell2mat(data(:,n));
    
    % Increase index
    s = t + 1;
end
    
% Set Inf values to NaN because
bpd(isinf(bpd)) = NaN;
    
% Set the current axes
axes(fig.ax); hold on;

% Jitter value
jitter = 0.4;

% Now draw a box plot
boxplot(bpd,{bpi(:,1) bpi(:,2)},...
    'Jitter',jitter,...
    'DataLim',[-Inf 1e6],...
    'LabelVerbosity','all',...
    'FactorSeparator',1);

set(gca,'FontSize',16);
set(findobj(gca,'Type','text'),...
    'FontSize',16,...
    'FontWeight','bold',...
    'VerticalAlignment','middle');


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%