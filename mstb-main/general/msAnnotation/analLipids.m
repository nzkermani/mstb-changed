function analLipids( op,ass,saveP )
%analLipids - a function for Nima's lipip analysis.
% James McKenzie, Imperial College, London, 2014.
%
%
% To do list:
%       x toggle between bar/box plots
%       x normalisation of peaks according to annotated species only
%       x sparse matrix denoting which files contain which histIDs
%       x export of 'processed raw' data
%       x relevant additions to the uiTree for this

% Decide if the savePath has been passed to the function, if not go with a
% new one...
if nargin == 0
    
    % Need to be able to reload lipid analyses...
    path = pwd;
    
    [name,path,~] = uigetfile({'*.mat'},'Select previous analLipids',path);
    if ~strcmp(path(end),filesep)
        path = [path filesep];
    end
    
    tmp = open([path name]);
    
    try   
        op  = tmp.nima.op;
        ass = tmp.nima.ass;
        sP  = [path name(1:end-4) '-' datestr(now,'yymmdd-HHMM') '.mat'];
    catch
        error('Cannot load this file');
    end
    
    
elseif nargin == 2
            
    % Determine file name...use the time of analysis
    fName = ['Lipids-' datestr(now,'yymmdd-HHMM') '.mat'];

    try
        sP = [op.paths.local 'MetaMerge' filesep fName];
    catch
        % Generate the path if none provided
        ppp = pwd;
        if ~strcmp(ppp(end),filesep)
            ppp = [ppp filesep 'LipidAnalysisResults' filesep];
        end

        % Create folder
        if ~exist(ppp,'dir')
            mkdir(ppp);
        end


        sP = [ppp fName];
        warning on all
        warning('Default save location');
    end    
elseif nargin == 3
    sP = saveP;    
end

% Add the input data to the guidata
op.XPeaks    = full(op.XPeaks);
op.XPeaksLog = log(full(op.XPeaks) + op.logOS);%full(op.XPeaksLog);
gui.op  = op;
gui.ass = ass;

% DO THESE FUNCTIONS JUST THE ONCE
% Start by drawing a generic window which will be filled with stuff when we
% decide what we want to include in it...
[gui.fig] = drawWindow;

% Now get the data for the table - this is the peak list essentially in a
% matlab table format... This needs to be done only once in reality.
[gui.tabPks] = tableData(op,ass);

% Here we will define the colour scheme to be used across the entire
% function
numG = numel(unique(gui.op.histID));
if numG == 2
    gui.cols = [1 0 0; 0 1 0];
elseif numG == 3
    gui.cols = [0 0 1; 1 0 0; 0 1 0];
else
    gui.cols = jet(numG);
end

% Actually put the data in the table
set(gui.fig.tabPks,...
    'Data',gui.tabPks.dat,...
    'ColumnName',gui.tabPks.colLabs,...
    'ColumnWidth',gui.tabPks.colWdth);

% Set the guidata for the window
guidata(gui.fig.fig,gui);

% Initialise the window with this catch-all function...
figPop(gui.fig.fig,[]);

% Now update the callback functions, which needs to be done only the once.
% Make it generic and provide no variables, as everything can be gathered
% from either guidata or the tables accessible via the handles
tbCallbacks(gui.fig.fig,[],sP)


% How many annotations?
numAnno = sum(cell2mat(gui.tabPks.dat(:,1)));
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' ...
    char(10) int2str(numAnno) ' annotations' char(10)...
    '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);

% Display save path location
disp(sP);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are the run only once functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = drawWindow
% Draw the main window

f0 = findobj('Name','Lipid Analysis');
close(f0);

fig.fig = figure('Name','Lipid Analysis',...
    'Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5],...
    'Menubar','none',...
    'Toolbar','figure',...
    'Number','off');

% Add a small spectral axes for the median spectra
axSz = [0.01 0.8 0.98 0.195];
fig.ax1 = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',axSz);

% Here a secondard axes for visualisation of results...
fig.ax2 = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.53 0.06 0.46 0.67],...
    'Visible','off');

% Three small tables down the LHS for the key information
fig.tabCls = uitable('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.01 0.51 0.48 0.23]);

fig.tabLng = uitable('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.01 0.26 0.48 0.23]);

fig.tabSat = uitable('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.01 0.01 0.48 0.23]);

fig.tabPks = uitable('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.51 0.01 0.48 0.73]);

% Buttons for the barcharts of the tables...
% fig.butCls = uicontrol('Parent',fig.fig,...
%     'Style','pushbutton',...
%     'Units','normalized',...
%     'Position',[0.46 0.51 0.03 0.05],...
%     'String','Bar');
% 
% fig.butLng = uicontrol('Parent',fig.fig,...
%     'Style','pushbutton',...
%     'Units','normalized',...
%     'Position',[0.46 0.26 0.03 0.05],...
%     'String','Bar');
% 
% fig.butSat = uicontrol('Parent',fig.fig,...
%     'Style','pushbutton',...
%     'Units','normalized',...
%     'Position',[0.46 0.01 0.03 0.05],...
%     'String','Bar');
    
[fig] = drawDBtoolbar(fig);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = drawDBtoolbar(fig)

% Let's add a new toolbar to it instead of the horrible push buttons
tb.tb = uitoolbar('Parent',fig.fig,...
    'Tag','lipidtb',...
    'Visible','on');

% Need to get the icons for the buttons...
icons = open(deSlash('msiUI/msiDB/misc/lipidIcons.mat'));

% Folder and save buttons
% tb.fold = uipushtool('Parent',tb.tb,...
%     'CData',icons.nima.fold,...
%     'ClickedCallback',{},...
%     'TooltipString','Default save folder',...
%     'Tag','uiEdit',...
%     'UserData',pwd);

tb.save = uipushtool('Parent',tb.tb,...
    'CData',icons.nima.save,...
    'ClickedCallback',{},...
    'TooltipString','Press me to save the results.',...
    'Tag','none',...
    'UserData',pwd);

tb.print = uipushtool('Parent',tb.tb,...
    'CData',icons.nima.print,...
    'ClickedCallback',{},...
    'TooltipString','Save the current figure...',...
    'Tag','none',...
    'UserData',pwd);

uipushtool('Parent',tb.tb,...
    'CData', NaN([50 50 3]),...
    'ClickedCallback',{},...
    'Separator', 'on');

% This is to show the table - depressed on initialisation
tb.tabPks = uitoggletool('Parent',tb.tb,...
    'CData',icons.nima.pks,...
    'State','On',...
    'ClickedCallback',{},...
    'TooltipString','Peak table view');

tb.barCls = uitoggletool('Parent',tb.tb,...
    'CData',icons.nima.cls,...
    'State','Off',...
    'ClickedCallback',{},...
    'TooltipString','Plot class',...
    'Tag','uiGraph',...
    'UserData','class');
tb.barLng = uitoggletool('Parent',tb.tb,...
    'CData',icons.nima.lng,...
    'State','Off',...
    'ClickedCallback',{},...
    'TooltipString','Plot length',...
    'Tag','uiGraph',...
    'UserData','length');
tb.barSat = uitoggletool('Parent',tb.tb,...
    'CData',icons.nima.sat,...
    'State','Off',...
    'ClickedCallback',{},...
    'TooltipString','Plot desaturations',...
    'Tag','uiGraph',...
    'UserData','desat');

uipushtool('Parent',tb.tb,...
    'CData', NaN([50 50 3]),...
    'ClickedCallback',{},...
    'Separator', 'on');

tb.grBox = uitoggletool('Parent',tb.tb,...
    'CData',icons.nima.box,...
    'State','On',...
    'ClickedCallback',{},...
    'TooltipString','Box plots',...
    'Tag','uiXXX',...
    'UserData','box');

tb.grBar = uitoggletool('Parent',tb.tb,...
    'CData',icons.nima.bar,...
    'State','Off',...
    'ClickedCallback',{},...
    'TooltipString','Bar charts',...
    'Tag','uiXXX',...
    'UserData','bar');

uipushtool('Parent',tb.tb,...
    'CData', NaN([50 50 3]),...
    'ClickedCallback',{},...
    'Separator', 'on');

tb.showall= uitoggletool('Parent',tb.tb,...
    'CData',icons.nima.ringfull,...
    'ClickedCallback',{},...
    'TooltipString','Full circle = all lengths',...
    'Tag','oddeven',...
    'State','on',...
    'UserData',false);

tb.showeven = uitoggletool('Parent',tb.tb,...
    'CData',icons.nima.ringdash,...
    'ClickedCallback',{},...
    'TooltipString','Dashed circle = even only',...
    'Tag','oddeven',...
    'State','off',...
    'UserData',true);

uipushtool('Parent',tb.tb,...
    'CData', NaN([50 50 3]),...
    'ClickedCallback',{},...
    'Separator', 'on');

tb.foldchange = uipushtool('Parent',tb.tb,...
    'CData',icons.nima.fold,...
    'ClickedCallback',{},...
    'TooltipString','Fold(er) changes',...
    'Tag','none',...
    'UserData',pwd);

uipushtool('Parent',tb.tb,...
    'CData', NaN([50 50 3]),...
    'ClickedCallback',{},...
    'Separator', 'on');

uipushtool('Parent',tb.tb,...
    'CData', NaN([50 50 3]),...
    'ClickedCallback',{},...
    'Separator', 'on');

tb.edit = uitoggletool('Parent',tb.tb,...
    'CData',icons.nima.locked,...
    'State','Off',...
    'ClickedCallback',{},...
    'TooltipString','When the key is shown the data is locked.',...
    'Tag','uiEdit',...
    'UserData',icons.nima.unlocked);





% What about enlarging the toolbar?
enlargeToolbar(tb.tb,50);

fig.tb = tb;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tab] = tableData(op,ass)
% Get the annotation data in a table for presentation...

% Quantity of mz values, hence table entries
numE = numel(op.cmz);

% An empty table
tab.dat = cell(numE,8);

% Cycle through the mz variables
for n = 1:numE
    
    
    % Can only focus on peaks that are singly annotated. This if statement
    % skips anything that is unknown.
    if isempty(strfind(ass.annoNam{n,1},'/')) && length(ass.annoNam{n,1}) > 3
        
        % This is the name of the annotated lipid
        tmp = ass.annoNam{n,1};
        
        % These are markers for finding the length/desat stuff
        cho = strfind(tmp,'(');
        chc = strfind(tmp,')');
        chl = strfind(tmp,':');

        % Need to exclude P-xx and O-xx lipids as I can't think of a good
        % way to display these in a grid. May also have to apply this to
        % Nima's analLipids function
%         qtmp = tmp(cho(1)+1:chl(1)-1);
%         if strcmp(qtmp(1:2),'P-') || strcmp(qtmp(1:2),'O-')
%             skip = true;
%         else
%             skip = false;
%         end

        if ~isempty(cho) && ~isempty(cho) && ~isempty(chl)
            
            % How many carbons?
            nCarb  = tmp(cho(1)+1:chl(1)-1);

            % This might include a P- or an O-...            
            if length(nCarb) > 2
                if strcmp(nCarb(1:2),'P-') || strcmp(nCarb(1:2),'O-')
                    nCarb = nCarb(3:end);
                end
            end
                        
            % Should be straight forward
            nDesat = tmp(chl(1)+1:chc(1)-1);
        end
        
        % To which class does this belong?
        annoCls = ass.annoCls{n,1};
    
        % What is the true mz of the assignment?
        fx  = find(ass.list(:,1) == n);
        ffx = ass.list(fx,2);         %#ok<FNDSB>
        trueMZ = ass.db.mz(ffx);
        
        
        % ppm difference
        ppmDev = 1e6 * ( op.cmz(1,n) - trueMZ) / trueMZ;
        ppmDev = sprintf('%0.2f',ppmDev);
        trueMZ = sprintf('%0.4f',trueMZ);
    else
        % Then there are multiple assignments, so do nothing
        nDesat = '';
        nCarb  = '';
        annoCls= '';     
        trueMZ = '';
        ppmDev = '';
    end
    
    % Now we need to look at the isotopic information
    if ass.isIso(n,1) == 0
        
        % Only display an M if we have identified an isotopic cluster
        nn = sum(nnz([ass.iso1(n,1) ass.iso2(n,2)]));
        if nn == 0
            isoStat = '';
        else
            isoStat = 'M';
                       
        end
        isoMZ = '';
        isoCorr = 0;
        isoRat = 0;
    else
        % Either M+0.5, M+1 M+2
        mzdiff = round((op.cmz(1,n) - ass.isIso(n,1)) / 0.5);
        if mzdiff == 1
            isoStat = 'M+1/2';
        elseif mzdiff == 2
            isoStat = 'M+1';
        elseif mzdiff == 4
            isoStat = 'M+2';
        else
            isoStat = 'ERR!';
            %mzdiff
            %op.cmz(1,n)
            %ass.isIso(n,1)
            %pause
        end
        
        % This is the mz value of the parent
        isoMZ = ass.isIso(n,2);
        
        % Calculate the correlation coefficient between this peak and
        % its parent M peak
        isoCorr = corr(op.XPeaks(:,n),op.XPeaks(:,ass.isIso(n,2)));
        %isoCorr = sprintf('%0.2f',isoCorr);
        
        isoA = op.XPeaks(:,n);
        isoB = op.XPeaks(:,ass.isIso(n,2));
        
        isoRat  = 100 * isoA ./ isoB;
        isoRat  = nanmean(isoRat);
        %isoStd  = nanstd(isoRat);
        %isoRat  = isoRat;
        
    end
    
    % This is the new complete logic for inclusion of all peaks, which
    % encompasses isotopic and assignment information. Start by default
    % with an excluded peak and look to find reasons to include it.
    tick = false;
    
    if ass.annoNum(n,1) == 1% && strcmp(ass.annoCls{n,1}(1:3),'GPL')
        % Thus only look at peaks with 1 assigned GLP candidate - it is a
        % little too complicated to manage multiple assignments.
        % By commenting the second condition in above, we now include the
        % free fatty acids in the calculations. If you don't want them,
        % then you'll need to limit the m/z range...
        
        switch isoStat
            
            case 'M'
                % This is included by default, just need to decide about
                % whether to include any of its so called isotopic peaks
                % but that is independent of this peak's inclusion
                tick = true;
                
            case 'M+1'
                
                % May be either the first or second isotope, dependent on
                % the the charge. Need to look into fixing this!
                if (op.cmz(1,n) - op.cmz(1,isoMZ)) > 0.75
                    % Only include if poor correlation to M
                    if isoCorr < 0.8 
                        tick = true;
                    elseif isempty(ass.annoNam{isoMZ,1})
                        tick = true;
                    end
                else
                    % Then this is the M+1 peak of an M+1/2 peak, thus
                    % treat the same as the M+2 peak
                    if isoRat > 20 && isoCorr < 0.5
                        tick = true;
                    end
                end
                
            case 'M+2'
                % Include if orrelation less than 0.5 and ratio greater
                % than 20%
                if isinf(isoRat)
                    
                elseif strcmp(isoCorr,'NaN')
                    
                else
                    
                    if isoRat > 20 && isoCorr < 0.5
                        tick = true;
                    end
                end                    
                
            otherwise
                % This is a peak sans isotopic friends. Should expect every
                % genuine peak to have a +1 isotope except those with
                % intensities below the limit of detection.
                tick = true;                
        end
        
    end

    % Convert isoRat into a string
    if isoRat == 0
        isoRat  = '';
        isoCorr = '';
    else
        isoRat  = [sprintf('%2.0f',isoRat ) '%'];
        isoCorr =  sprintf('%0.2f',isoCorr);
    end
    
    
    
    
    % Add the data here once you've done the necessary calculations    
    tab.dat{n,1} = tick;
    tab.dat{n,2} = op.cmz(1,n);
    tab.dat{n,3} = isoStat;
    tab.dat{n,4} = isoMZ;
    tab.dat{n,5} = isoCorr;
    tab.dat{n,6} = isoRat;
    tab.dat{n,7} = ass.annoNum(n,1);
    tab.dat{n,8} = ass.annoNam{n,1};
    tab.dat{n,9} = trueMZ; %ass.annoMZ(1,n);
    tab.dat{n,10} = ppmDev;
    tab.dat{n,11} = annoCls;
    tab.dat{n,12} = nCarb;
    tab.dat{n,13} = nDesat;
        
end

tab.colLabs = {'?','m/z','Iso','M','Corr','Rat','nID','ID','True m/z','±ppm','Class','nCarb','nDesat'};
tab.colWdth = { 20    60    35  30     40    40    25  100         60     50      80      50       50}; 

% Finally calculate the spectral sum of annotated peaks! This is best as a
% separate function at is will change when the table is edited to
% include/exclude peaks...
%[tab.int] = calcSumAnnotations(op,tab.dat);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figPop(src,event)
% Main function for the execution of the functions within this ui

gui = guidata(src);

% Calculate the normalised spectral intensities
calcSumAnnotations(src,[]);

% Let's add in the median spectra
addNormedSpectra(src,[]);

% Update the tables containing class/length/desat
annoInfo(src,[])

% Add the triangle markers on the top axes
addAnnotations(src,[]);

% Statistical tests
testStats(src,[],gui.fig.tabCls,'cls');
testStats(src,[],gui.fig.tabLng,'lng');
testStats(src,[],gui.fig.tabSat,'sat');

% So now it would appear that all of the data acquisition functions and
% calculations have been performed. Now need to return to the main function
% and update callbacks on buttons and stuff

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tbCallbacks(src,~,sP)
% Update all of the callbacks for the buttons in the toolbar

gui = guidata(src);

% This is the table toggle button
set(gui.fig.tb.tabPks,'ClickedCallback', {@showPeakTable,gui.fig});

% These 'toggle' betwen the graphs
set(gui.fig.tb.grBox,'ClickedCallback',{@setGraphChoice,gui.fig.tb.grBar});
set(gui.fig.tb.grBar,'ClickedCallback',{@setGraphChoice,gui.fig.tb.grBox});

% Logic for the graph buttons
set(gui.fig.tb.barCls,'ClickedCallback', {@plotButtonLogic,gui.fig});
set(gui.fig.tb.barLng,'ClickedCallback', {@plotButtonLogic,gui.fig});
set(gui.fig.tb.barSat,'ClickedCallback', {@plotButtonLogic,gui.fig});

% Table edit function
set(gui.fig.tb.edit,'ClickedCallback', {@tableEdit});%,fig,medSpec,class,op.histID});

% Save function 
set(gui.fig.tb.save,'ClickedCallback', {@saveResults,sP});

% Print figure function
set(gui.fig.tb.print,'ClickedCallback',{@printFigure,sP});

%set(gui.fig.tb.showall, 'ClickedCallback', {@oddEvenLogic,gui.fig});
%set(gui.fig.tb.showeven,'ClickedCallback', {@oddEvenLogic,gui.fig});

set(gui.fig.tb.showall, 'ClickedCallback', {@setGraphChoice,gui.fig.tb.showeven});
set(gui.fig.tb.showeven,'ClickedCallback', {@setGraphChoice,gui.fig.tb.showall});

% For the fold changes...
set(gui.fig.tb.foldchange,'ClickedCallback', {@foldChange,sP});

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are the run many times functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calcSumAnnotations(src,event)
% Sum up the intensities of the annotated peaks in each file and return
% these...

gui = guidata(src);

sumInt = zeros(size(gui.op.XPeaks,1),1);
szTab = size(gui.tabPks.dat,1);
for n = 1:szTab
    tmp = gui.tabPks.dat(n,1);
    if tmp{:} %#ok<BDSCA>
        sumInt(:,1) = sumInt(:,1) + gui.op.XPeaks(:,n);        
    end    
end

% Update guidata and return
gui.sumInt = sumInt;
guidata(src,gui);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function addNormedSpectra(src,event)%fig,op,specSum)
% Calculate the median spectra from op, and add it to the graph...

normMethod = 'none'; % TIC LIC none

% Gather guidata
gui = guidata(src);

% Here is where we normalise the XPeaks to the sum of the annotated
% species. This is an improvement on the earlier version which normalised
% based on all species within the range.
switch lower(normMethod)    
    case 'lic'
        meds = 1000 * bsxfun(@rdivide,gui.op.XPeaks,gui.sumInt);        
    case 'tic'
        meds = 1000 * bsxfun(@rdivide,gui.op.XPeaks,sum(gui.op.XPeaks,2));        
    case 'none'
        meds = gui.op.XPeaks;
    otherwise
        error('there is no otherwise');
end

% Determine the unique histIDs
[unq,~,~] = unique(gui.op.histID);
    
% Formulate the colour scheme
cols = gui.cols;

% Add to the axes...    
axes(gui.fig.ax1); 
cla reset;
hold on;

% Insert zeros into the spectral for the plot to make it look nicer
[new.mz,new.sp] = insertZeros(gui.op.cmz,meds);
drawMethod = 'plot';

h = zeros(numel(unq),1);
for n = 1:numel(unq)
    
    fx = strcmp(gui.op.histID,unq{n});
    fx = find(fx == 1);
    
    switch drawMethod
        case 'stem'
            for r = 1:numel(fx)
                %tmp = plot(op.cmz,meds(fx(r),:),'Color',cols(n,:));
                tmp = stem(gui.op.cmz,meds(fx(r),:),...
                    'LineWidth',1,...
                    'Color','k',...%cols(n,:),...
                    'MarkerFaceColor',cols(n,:),...
                    'MarkerEdgeColor','k');%cols(n,:));

                if r == 1
                    h(n,1) = tmp(1);
                end
            end
            
        case 'plot'
            tmp = plot(new.mz,new.sp(fx,:),...
                'LineWidth',1,...
                'Color',cols(n,:));
            h(n,1) = tmp(1);
    end
    
end
leg = legend(h,unq);
set(leg,'Color','none');

box on;
set(gca,...%'Color','none',...
    'YTick',[]);

% Update guidata and return
gui.spec = meds;
gui.unqHist = unq;
guidata(src,gui);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function annoInfo(src,event)
% This function takes the information in tabPks and produces the individual
% tables containing information about class/length/desat

% Get all the pertinent guidata
gui = guidata(src);

% Get the table data...
dat = gui.tabPks.dat;       %get(fig.tabPks,'Data');
cnm = gui.tabPks.colLabs;   %get(fig.tabPks,'ColumnName');

% Get the normalised spectra medSpec
medSpec = gui.spec;

% Start with the length...
lngind = find(strcmp(cnm,'nCarb') == 1);
lng = str2double(dat(:,lngind));
lng(isnan(lng)) = 0;
unqLng = unique(lng);
unqLng = unqLng(2:end)'; % remove zero...

% Now do the ndesat
satind = find(strcmp(cnm,'nDesat') == 1);
sat = str2double(dat(:,satind));
sat(isnan(sat)) = 0;
unqSat = unique(sat)';

% Now the class of lipids - this contains everything, but zero-valued
% classes can be removed from the table later on.
clsind = find(strcmp(cnm,'Class') == 1);
cls = dat(:,clsind);
unqCls = unique(cls);
unqCls = unqCls(2:end)';

% Tables for the data
tabLng = zeros(size(medSpec,1),numel(unqLng));
tabSat = zeros(size(medSpec,1),numel(unqSat));
tabCls = zeros(size(medSpec,1),numel(unqCls));

for n = 1:size(dat,1)
    
    if cell2mat(dat(n,1)) % true, then can include...
        
        % What are its classwise median intensities?
        val = medSpec(:,n);
        
        % What is the length of this annotation?
        tmpLng = str2double(dat(n,lngind));
        
        % To which length 'bin' does this correspond?
        lngBin = find(unqLng == tmpLng);
        
        % Add into the table / thing
        tabLng(:,lngBin) = bsxfun(@plus,tabLng(:,lngBin),val);
        
        
        % Now do the same for the number of desaturations...
        tmpSat = str2double(dat(n,satind));
        satBin = find(unqSat == tmpSat);
        tabSat(:,satBin) = bsxfun(@plus,tabSat(:,satBin),val);
        
        % Now finally let's do the same for the class of the lipids,
        % hopefully not too difficult...
        tmpCls = dat(n,clsind);
        clsBin = strcmp(unqCls,tmpCls);
        clsBin = find(clsBin == 1);
        tabCls(:,clsBin) = bsxfun(@plus,tabCls(:,clsBin),val);
    end
end

% Ditch the extra entries from the tables, i.e. no zeros
[~,fx] = find(sum(tabCls,1) == 0);
tabCls(:,fx) = [];
unqCls(:,fx) = [];

[~,fx] = find(sum(tabSat,1) == 0);
tabSat(:,fx) = [];
unqSat(:,fx) = [];

[~,fx] = find(sum(tabLng,1) == 0);
tabLng(:,fx) = [];
unqLng(:,fx) = [];


[b,i] = sort(gui.op.histID);

tabCls = tabCls(i,:);
tabLng = tabLng(i,:);
tabSat = tabSat(i,:);

histID = b;
try
    newPID = gui.op.pID(:);
catch
    newPID = gui.op.patientID(:);
end

% Format into tabular format...
set(gui.fig.tabLng,'Data',tabLng,...
    'ColumnName',unqLng,...
    'RowName',histID);

set(gui.fig.tabSat,'Data',tabSat,...
    'ColumnName',unqSat,...
    'RowName',histID);

set(gui.fig.tabCls,'Data',tabCls,...
    'ColumnName',unqCls,...
    'RowName',histID);

% We won't reutrn the tables to the guidata as there is no need,
% considering that they can be accessed via the table handle
gui.labels.cls = unqCls;
gui.labels.lng = unqLng;
gui.labels.sat = unqSat;
gui.labels.pID = newPID;
guidata(src,gui);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function addAnnotations(src,event)
% Add annotations to the top most plot, showing subclasses so see how much
% grouping actually goes on.

% Gather the guidata
gui = guidata(src);

% Base this on the ticked entries in the table...
tab = gui.tabPks.dat;
cnm = gui.tabPks.colLabs;

clsind = find(strcmp(cnm,'Class') == 1);

% The numbers of classes that are annotated!
labs = get(gui.fig.tabCls,'ColumnName');
numL = numel(labs);

% Define the colours
cols = hsv(numL);

% Get the ylim so we know how high to plot it...
yl = ylim;
yl = yl(2);

xl = xlim;

% This is for the elevation
split = 0.4 / numL;
yup = 0.55:split:0.95;

% How wide to make each band?
marg = 4;

for n = 1:size(tab,1)
    
    if cell2mat(tab(n,1))
        
        mz = cell2mat(tab(n,2));
        cl = tab(n,clsind); %#ok<FNDSB>
        
        ind = strcmp(labs,char(cl));
        
        if sum(ind) == 0
            continue;
        end
        
        fx = find(ind == 1);
  
        ypos = [yup(fx) yup(fx+1) yup(fx+1)] * yl;        
        patch([mz mz+marg mz-marg],ypos,...
            cols(fx,:),'FaceAlpha',0.125,...
            'EdgeColor','none');
        
    end    
end

% Text labels down the LHS to tell us which colour patch is for which
% class
for n = 1:numL
    text(xl(1)+10,mean(yup(n:n+1))*yl,labs{n},'VerticalAlignment','middle')
end

xlim(xl);

% Don't need to update guidata as nothing was changed

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testStats(src,~,fig,name)
% Calculate the test statistics for the groups...

% Gather the appropriate data...
gui = guidata(src);

dat = get(fig,'Data');
row = get(fig,'RowName');
col = get(fig,'ColumnName');

% How many distinct groups?
[grps,~,inds] = unique(row);
numG          = numel(grps);

% How many variables?
numV = size(dat,2);

% If there are two groups, do MannWhitneyU test
if numG == 2
    
    fx = find(inds == 1);
    fy = find(inds == 2);
    
    mwu = zeros(numV,2);
    for n = 1:numV
        [mwu(n,1),mwu(n,2),~] = ranksum(dat(fx,n),dat(fy,n)); %#ok<FNDSB>
    end
    
else
    mwu = zeros(numV,2);    
end

% Now do anova in cases of there being more than one group (note that post
% hoc analysis is not performed)
ano = zeros(numV,1);
pha = ones(numG,numG,numV);
for n = 1:numV
    [ano(n,1),~,stts] = anova1(dat(:,n),row,'off');
    
    % Posthoc ANOVA if ano(n,1) is significant
    if ano(n,1) <= 0.05
        
        mc = multcompare(stts,'ctype','hsd','display','off');        
        for r = 1:size(mc,1)
            pha(mc(r,1),mc(r,2),n) = mc(r,6);
            pha(mc(r,2),mc(r,1),n) = mc(r,6);
        end        
    end
end

% Kruskal Wallis testing here
krk = zeros(numV,1);    
phk = ones(numG,numG,numV); % posthoc KW if significant
for n = 1:numV
    [krk(n,1),~,~] = kruskalwallis(dat(:,n),row,'off');
    
    % If the KW was significant AND there are more than two groups, do
    % individual MWU tests between each group. You can worry about multiple
    % sampling later on!
    if krk(n,1) <= 0.05 && numG > 2        
        for j = 1:numG
            for k = 1:numG
                fx = find(inds == j);
                fy = find(inds == k);                
                phk(j,k,n) = ranksum(dat(fx,n),dat(fy,n));
            end
        end
    else
        % do nothing
    end
        
end

varNam = {'MannWhitneyU','ANOVA','KruskalWallis'};

% Display in a table...
results = table(mwu(:,1),ano(:,1),krk(:,1),...
    'RowNames',cellstr(col),'VariableNames',varNam) %#ok<NOPRT>

% Update the guidata with the results
gui.stats.(name)    = results;
gui.statph.(name)   = phk;
gui.statan.(name)   = pha;

guidata(src,gui);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showPeakTable(src,~,fig)
% Just show the peak table when this is turned on, nothing when turned off

% Get state of this button
state = get(src,'State');

switch state
    case 'on'
        
        % Make the axes invisible
        set(fig.ax2,'Visible','off');
        
        % Make the table visible
        set(fig.tabPks,'Visible','on');
        
        % Find the handles of the 3 graph buttons and turn off
        f0 = findobj('Tag','uiGraph');
        set(f0,'State','off');        
        
    case 'off'        
        
        % This button cannot be deselected itself
        set(src,'State','on');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotButtonLogic(src,~,fig)
% Decide which buttons need to be deselected and which type of plot to
% draw, all too complicated!

% First decide if this button is already pressed
state = get(src,'State');

if strcmp(state,'on')
    
    % Disable all other graph push buttons...
    f0 = findobj('Tag','uiGraph');
    set(f0,'State','off');
    
    % And the table peak button
    set(fig.tb.tabPks,'State','off');
    
    % Re-enable the clicked button
    set(src,'State','on');

    % Get the choice of box or
    if strcmp(get(fig.tb.grBox,'State'),'on')
        plotMethod = 'box';
    else
        plotMethod = 'bar';
    end
    
    % Odd/even logic?
    f1 = findobj('Tag','oddeven','State','on');
    oeType = get(f1,'UserData');
    
    % Now call the plotting function...
    plotType = get(src,'UserData');
    
    doChartPlot(src,[],plotType,plotMethod,oeType);

else    
    % Show the table...
    set(fig.ax2,'Visible','off');
    set(fig.tabPks,'Visible','on');
    set(fig.tb.tabPks,'State','on');
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGraphChoice(src,~,src2)
% Set the bar/box chart choice

state = get(src,'State');
switch state    
    case 'on'        
        set(src2,'State','Off');        
    case 'off'        
        set(src2,'State','on');        
end

% Also can call the graph function here, assuming that one of the plot
% buttons is toggled
f0 = findobj('Tag','uiGraph','State','on');
f1 = findobj('Tag','oddeven','State','on');

% This is the plot type, i.e. from which table is the data coming?
if ~isempty(f0)
    plotType = get(f0,'UserData');
    
    % odd/even?
    oeType = get(f1,'UserData');
    
    f0 = findobj('Tag','uiXXX','State','on');
    plotMethod = get(f0,'UserData');
    
    % Call the plot function here
    doChartPlot(src,[],plotType,plotMethod,oeType);
    
else
    % Do nothing, as the table is shown
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doChartPlot(src,~,style,type,oddEven)%,fig,hTab,~,tstTab)
% This is the plotting function for displaying either the box/bar chart for
% each of the three tables 

% Gather the guidata
gui = guidata(src);

% Decide which of the data table's we are plotting
switch style
    case 'class'
        hTab = gui.fig.tabCls;
        tstTab = gui.stats.cls;
        xLab = 'Lipid class';
    case 'length'
        hTab = gui.fig.tabLng;
        tstTab = gui.stats.lng;
        xLab = 'Lipid length';
    case 'desat'
        hTab = gui.fig.tabSat;
        tstTab = gui.stats.sat;
        xLab = 'Lipid desaturations';
end

% Are we going to plot only the even values?
if strcmp(xLab,'Lipid length') && oddEven %isTrue
    special = true;
else
    special = false;
end


% Get the data from the table's handle
td = get(hTab,'Data');
cn = cellstr(get(hTab,'ColumnName'));
rn = cellstr(get(hTab,'RowName'));

% Here we do a little trimming
if special
    % Convert the cell to double
    abc = str2num(cell2mat(cn));
    
    % Keep indices of even values
    keep = ~mod(abc,2) & abc >= 32;
    
    % Trim the tables...
    td = td(:,keep);
    cn = cn(keep);
    
end

% Define the colours...
unq = unique(rn);
cols = gui.cols;

% Where to plot this thing?
set(gui.fig.tabPks,'Visible','off');
set(gui.fig.ax2,   'Visible','on');

% Reset and hold 
axes(gui.fig.ax2);
cla reset
hold on;

% Need to arrange the data into a single longer matrix, i.e. cancer in the
% first set of columns, then healthy concatenated afterwards
numO = zeros(numel(unq),1);
for n = 1:numel(unq)    
    fx = strcmp(rn,unq{n});
    fx = find(fx == 1);
    numO(n,1) = numel(fx);    
end

% How many columns are there to be? 
unqCols = unique(cn);
numC = numel(unqCols) * numel(unq);

% Empty matrix of NaNs...box plot data
bpd = NaN(max(numO),numC);

% Now place into the NaN matrix
colM = [1 numel(unqCols)];
for n = 1:numel(unq)    
    % Find observations for this group, unq{n}
    fx = strcmp(rn,unq{n});
    fx = find(fx == 1);    
    bpd(1:numO(n,1),colM(1):colM(2)) = td(fx,:);     %#ok<FNDSB>
    colM = colM + numel(unqCols);
end



% We've previously done the anova/mu tests, so let's just decide which are
% significant in advance
symb = {'*','**','***'};
symu = {'u','uu','uuu'};
symk = {'k','kk','kkk'};
sigv = [0.05 0.01 0.001];

% Significant values from ANOVA
ano = bsxfun(@le,tstTab.ANOVA,sigv);
ano = sum(ano,2);

krk = bsxfun(@le,tstTab.KruskalWallis,sigv);
krk = sum(krk,2);

% If we've done MWU need to decide which are more significant than others
if numel(unq) == 2
    mwu = bsxfun(@le,tstTab.MannWhitneyU,sigv);
    mvw = sum(mwu,2);
end

% May need to trim out the significant values from the stat tests
if special
    ano = ano(keep,:);
    krk = krk(keep);
    if numel(unq) == 2
        mwu = mwu(keep);
        mvw = mvw(keep);        
    end
end

% Here loop between either of the choices of plots
switch type    
    case 'bar'
        
        % Now let's calculate the mean of each group for the bar chart
        grpmean = nanmean(bpd,1);
        
        % And also the standard error
        stdev  = nanstd(bpd,[],1);
        %no     = sum(~isnan(bpd),1);
        %stderr = stdev ./ sqrt(no);

        % This reshapes the data to make it one column per tissue type,
        % e.g. Cancer(:,1) / Healthy(:,2)
        grpmean = reshape(grpmean,numel(grpmean)/numel(unq),numel(unq));
        stdev   = reshape(stdev,  numel(grpmean)/numel(unq),numel(unq));
        %stderr  = reshape(stderr, numel(grpmean)/numel(unq),numel(unq));
        
        % Now lets' do the bar chart drawning
        h = bar(grpmean,'grouped');

        % Now change colours and add the error bars...
        for n = 1:numel(unq)
            set(h(n),'FaceColor',cols(n,:),'EdgeColor','none');
            
            % X location of the error bars?
            xpos = mean(get(get(h(n),'Children'),'XData'),1)';
            ypos = grpmean(:,n);
            
            zpos1 = grpmean(:,n) + stdev(:,n);
            %zpos2 = grpmean(:,n) + stderr(:,n);

            % Draw the lines on the plot...            
            la = [xpos xpos];
            lb = [ypos zpos1];
            line(la',lb','LineWidth',1,'Color',cols(n,:));            
        end
        

        % Now add in the ANOVA / MWU markers of significance
        
        % These are the x poistions...
        xpos = 1:numel(unqCols);
        ypos = ylim; 
        ypos = ypos(2);
        
        % Now add the symbols where necessary...
        for n = 1:numel(ano)            
            
            % Add in anova
            if ano(n,1) > 0                
                text(xpos(n),ypos*0.9,symb{ano(n,1)},'FontSize',20,...
                    'VerticalAlignment','middle',...
                    'HorizontalAlignment','center',...
                    'ButtonDownFcn',{@showPostHoc,n,style,special,cn,'an'});
            end
            
            % Add in MWU
            if numel(unq) == 2 && mvw(n,1) > 0
                text(xpos(n),ypos*0.95,symu{mvw(n,1)},'FontSize',12,...
                    'VerticalAlignment','middle',...
                    'HorizontalAlignment','center');
            end
            
            % Add in Kruskal
            if krk(n,1) > 0
                text(xpos(n),ypos*0.85,symk{krk(n,1)},'FontSize',12,...
                    'VerticalAlignment','middle',...
                    'HorizontalAlignment','center',...
                    'ButtonDownFcn',{@showPostHoc,n,style,special,cn,'kw'});
            end

        end
        
        % Set the x axes limit nicely
        xlim([0.5 size(cn,1)+0.5]);

        % This is font size and grouping informatino
        set(gca,'FontSize',10,...
            'XTick',1:numel(cn),...
            'XTickLabel',cn);
        
        ylabel('Mean intensity (with standard deviation bar)');


    case 'box'
        
        % Calculate grouping variable
        c1 = repmat(cn,numel(unq),1);
        c2 = cell(numel(c1),1);
        colM = [1 numel(unqCols)];
        for n = 1:numel(unq)
            tmp = repmat({unq{n}},numel(unqCols),1); %#ok<CCAT1>
            c2(colM(1):colM(2)) = tmp;
            colM = colM + numel(unqCols);
        end
        cc = {c1 c2};
        
        % Draw the boxplot
        boxplot(bpd,cc,...
            'factorseparator',1,...
            'factorgap',[1 0.5],...
            'labelverbosity','minor',...
            'labelorientation','inline',...
            'colors',cols,...
            'boxstyle','filled');
        

        % Add in the ANOVA / MWU stat markers...
        
        % These are the x poistions...
        xvals = xlim;
        numBP = size(bpd,2) / numel(unq);
        sep = (xvals(2) - xvals(1)) / (numBP+1);
        xpos = xvals(1):sep:xvals(2);
        xpos = xpos(2:end-1);%- sep/2;
        
        
        %xpos = 1:numel(unqCols) / numel(unq);
        ypos = ylim; 
        ypos = ypos(2);
        
        % Now add the symbols where necessary...
        for n = 1:numel(ano)            
            if ano(n,1) > 0                
                text(xpos(n),ypos*0.9,symb{ano(n,1)},'FontSize',20,...
                    'VerticalAlignment','middle',...
                    'HorizontalAlignment','center',...
                    'ButtonDownFcn',{@showPostHoc,n,style,special,bpd,'an'});
            end
            
            if numel(unq) == 2 && mvw(n,1) > 0
                text(xpos(n),ypos*0.95,symu{mvw(n,1)},'FontSize',12,...
                    'VerticalAlignment','middle',...
                    'HorizontalAlignment','center');
            end
            
            if krk(n,1) > 0
                text(xpos(n),ypos*0.85,symk{krk(n,1)},'FontSize',12,...
                    'VerticalAlignment','middle',...
                    'HorizontalAlignment','center',...
                    'ButtonDownFcn',{@showPostHoc,n,style,special,bpd,'kw'});
            end
            
        end        
    otherwise
        disp('there is no otherwise');
end

% Plonk in an x axis label
xlabel(xLab);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tableEdit(src,~)%,fig,medSpec,class,histID)
% Allow the table to be edited, perhaps enlarging it too...

% Gather the guidata
gui = guidata(src);

% Get the state and decide if we are turning on / off
state = get(src,'State');

switch state
    case 'on'
        
        % When the button is turned on, then the table becomes active for
        % editing, but not all of it should be editable.
        
        % First need to display the table and ditch the graphs...
        set(gui.fig.ax2,   'Visible','off');
        set(gui.fig.tabPks,'Visible','on');
        
        % Disable graph buttons
        f0 = findobj('Tag','uiGraph');
        set(f0,'State','off','Enable','off');        
        set(gui.fig.tb.tabPks,'State','on');
               
        % Change the key to the pencil for aesthetic effect
        tmp = get(src,'UserData');
        set(src,'Userdata',get(src,'CData'));
        set(src,'CData',tmp);
        
        % Set the table to be editable
        set(gui.fig.tabPks,'ColumnEdit',...
            [true false false false false ...
            false true true false false ...
            true true true]);
        
    case 'off'
        
        
        % Change image back to the key...
        tmp = get(src,'UserData');
        set(src,'Userdata',get(src,'CData'));
        set(src,'CData',tmp);

        % Re-enable the push buttons
        f0 = findobj('Tag','uiGraph');
        set(f0,'State','off','Enable','on');        

        % Change table editabilitiy to fully false
        set(gui.fig.tabPks,'ColumnEdit',false);      
        
        % Add the new table data to the structure
        newDat = get(gui.fig.tabPks,'Data');
        gui.tabPks.dat = newDat;
        
        set(gui.fig.tabPks,'Data',newDat);
        
        % Now run the figPop function which should repopulate the tables of
        % the figure with the new information
        guidata(src,gui);
        
        figPop(src,[]);
        
end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveResults(src,~,sP)
% Save the required results to the workspace
% assignin('base','desiredName',CurrentName);

% Gather the guidata
gui = guidata(src);

% Here a sparse matrix showing which histIDs are in which sample

% Find unique sampleIDs
unqSamp = unique(gui.op.sampleID);

% How many unique histIDs?
unqHist = unique(gui.op.histID);

sp = cell(numel(unqSamp),numel(unqHist));

for n = 1:numel(unqSamp)
    
    % These are the entries for this sample (perhaps more than one histID)
    fx = strcmp(gui.op.sampleID,unqSamp{n});
    fx = find(fx == 1);
    
    for r = 1:numel(unqHist)
        sp{n,r} = 0;
    end
    
    for r = 1:numel(fx)        
        % This is the actual name of the histID
        hid = gui.op.histID{fx(r)};
        
        % Convert to the number based system
        hid = strcmp(unqHist,hid);
        hid = find(hid == 1);
        sp{n,hid} = 1; %#ok<FNDSB>
    end
end

tab = cell2table(sp,'VariableNames',unqHist,'RowNames',unqSamp);
tab2 = ['File' unqHist'];
tab2 = cat(1,tab2,[unqSamp sp]);

% Return it to the workspace
assignin('base','laClass', tab2);

% Now export the table data for class/length/desat and the labels for rows
% and columns. Then write another function for the boxplotting according to
% the metadata...
mrg.labs = gui.labels;
mrg.cls  = get(gui.fig.tabCls,'Data');
mrg.lng  = get(gui.fig.tabLng,'Data');
mrg.sat  = get(gui.fig.tabSat,'Data');
assignin('base','laMRG',mrg);


% Now export the raw processed data with appropriate information? Just the
% annotated variables essentially. We will export the 'raw' data and the 
% sum normalised data made within this function.
vars.histID = gui.op.histID;
vars.sampID = gui.op.sampleID;


% Need to find out which variables are annotated? Best to do this from the
% tabPks
tab = get(gui.fig.tabPks,'Data');

% Get the ticked ones
ticked = cell2mat(tab(:,1));
[fx,~] = find(ticked == 1);

% Now get the actual variables
vars.raw = gui.op.XPeaks(:,fx);
vars.pro = gui.spec(:,fx);

vars.mz  = gui.op.cmz(1,fx);

% Need the m/z and annotation for these variables
%vars.mz = tab(fx,2);
vars.anno = tab(fx,8);
vars.cls = tab(fx,11);
vars.lng = str2num(char(tab(fx,12))); %#ok<*ST2NM>
vars.sat = str2num(char(tab(fx,13)));

% Statistical test tables
vars.stats.cls = gui.stats.cls;
vars.stats.lng = gui.stats.lng;
vars.stats.sat = gui.stats.sat;

% Return to the workspace for analysis by Nima
assignin('base','laVars',vars);

% Define a nice picture
cdata = rand(10,10,3);
map = hsv(10);

% Save the results...
laVars = vars; %#ok<*NASGU>
laClass= tab;
save(sP,'laVars','laClass')

% Write a lovely message
msgbox({'Results returned to the workspace:' 'laClass' 'laVars' sP},...
    'Save Success',...
    'custom',cdata,map,'modal');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printFigure(src,~,sP)
% Save the figure to a file, having opened it in a separate window just for
% plotting...

% Get the guidata
gui = guidata(src);

% Prepare default path name
slsh = strfind(sP,filesep);
defPath = sP(1:slsh(end));

% Ask the user for a file name...
[fN,fP,fI] = uiputfile(...
    {'*.pdf';'*.png';'*.jpg'},...
    'Save As...',defPath);

% Progress for wait bar
wb = waitbar(1/2,'Printing figure');

% Various options...
switch fI    
    case 1 % pdf        
        export_fig(gui.fig.ax2,[fP fN], '-m2', '-transparent');
        
    otherwise
        export_fig(gui.fig.ax2,[fP fN], '-m2');
end

% Delete waitbar
delete(wb);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showPostHoc(src,event,i,style,specialTrim,labs,uniMethod)

% Here gather the guidata
gui = guidata(src);

% Decide which of the graphs is showing...
switch uniMethod
    case 'kw'
        umLab = 'Kruskal-Wallis test with Mann-Whitney posthoc';   
        
        switch style
            case 'class'
                tab = gui.statph.cls;
                lab = gui.labels.cls{i};
            case 'length'
                tab = gui.statph.lng;   
                lab = int2str(gui.labels.lng(i));
            case 'desat'
                tab = gui.statph.sat;
                lab = int2str(gui.labels.sat(i));
        end
        addMess = 'NO regard has been given to the repeated measures problem.';

        
    case 'an'
        umLab = 'ANOVA test with Tukey''s HSD posthoc';      

        switch style
            case 'class'
                tab = gui.statan.cls;
                lab = gui.labels.cls{i};
            case 'length'
                tab = gui.statan.lng;   
                lab = int2str(gui.labels.lng(i));
            case 'desat'
                tab = gui.statan.sat;
                lab = int2str(gui.labels.sat(i));
        end
        addMess = 'Tukey''s HSD at \alpha = 0.05';
        
    otherwise
        % there is no otherwise
end

% Need to get the status of the trimmed length box as this will disrupt the
% table as the rows won't match...
if specialTrim
    warning('Stats not shown when in EVEN only mode');
    return
end

% Put into table
try
    t = array2table(tab(:,:,i),'RowNames',gui.unqHist,'VariableNames',gui.unqHist);
catch
    t = array2table(tab(:,:,i),'RowNames',gui.unqHist);
end

% Make a nice output showing everything that is pertinent
disp('Posthoc Analysis Tests...');
disp(['... ' umLab ' ...']);
disp(['...for ' lab]);
disp(t)
disp(addMess);
disp(['*****************' char(10) char(10)])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function foldChange(src,event,sP)
% Figures for the displayment of fold changes

% Get the guidata
gui = guidata(src);

f0 = findobj('Tag','meanmedianfoldchange');
close(f0);

% Indices and labels for the histIDs
[hID,~,hInd] = unique(gui.op.histID);

% Ask the user which groups to fold
[g1,g2,lng,qv] = foldChangeDialog(hID);

indX = hInd == g1;
indY = hInd == g2;
hID = {hID{g1}; hID{g2}};

% Run through each entry in the table and put its fold change (to be
% calculated) into the appropriate slot
numT = size(gui.tabPks.dat,1);
fcInfo =  cell(numT,7);
fcVals = zeros(numT,3); % mean, median, MWU-p, MWU-q
mwu = zeros(numT,2);    % MWU p and q values
ano = zeros(numT,2);

for n = 1:numT
    
    % Determine its three vital stats...
    tC = gui.tabPks.dat{n,11};
    tL = str2num(gui.tabPks.dat{n,12});
    tS = str2num(gui.tabPks.dat{n,13});
        
    % Only continue if it is annotated
    if ~isempty(tC) && tL >= min(lng) && tL <= max(lng) %&& strcmp(tC(1:3),'GPL')
                
        % Calculate the fold changes
        meanFC   =   mean(gui.spec(indX,n)) /   mean(gui.spec(indY,n));
        medianFC = median(gui.spec(indX,n)) / median(gui.spec(indY,n));
        
        % Calculate the MWU value for the difference betweeen groups
        mwu(n,1) = ranksum(gui.spec(indX,n),gui.spec(indY,n));    
        
        % Do ANOVA perhaps instead of / as well as t-testing
        ano(n,1) = anovan(gui.spec(:,n),hInd,'display','off');
        
        
        % Parse out the GPL stuff before the comma
        fx = strfind(tC,',');
        
        % Determine if it is a plasmalogen?
        plP = strfind(gui.tabPks.dat{n,8},'P-');
        plO = strfind(gui.tabPks.dat{n,8},'O-');
        
        if isempty(plP) && isempty(plO)
            isPlas = '';
        elseif isempty(plP)
            isPlas = 'O';
        elseif isempty(plO)
            isPlas = 'P';
        else
            isPlas = 'PO!';
            warning('Something seriously wrong here!');
        end
        
        % Plonk into the structure / cell
        fcInfo{n,1} = tC(fx(1)+1:end);
        fcInfo{n,2} = tL;
        fcInfo{n,3} = tS;
        fcInfo{n,4} = log2( meanFC );
        fcInfo{n,5} = log2(medianFC);  
        fcInfo{n,6} = n;
        fcInfo{n,7} = isPlas;
        
        %fcVals(n,1) = log2( meanFC );
        %fcVals(n,2) = log2(medianFC);        
    end
end


% Calculate the p-values for the t-test (MWU values done previously)
%ttp = zeros(numT,2);
%ttp(:,1) = mattest(gui.spec(indX,:)',gui.spec(indY,:)');
% Now we are using one-way ANOVA to calucalte the values, rather than the
% t-testing, but the variable shall remain named as ttp
ttp = ano;

% How do we determine the q-values for the statistical tests?
switch qv{1}
    case 'bhy'
        % The Benjamini–Hochberg–Yekutieli method, apparently...
        [~,ttp(:,2)] = getBHYqVls(ttp(:,1)',qv{2});
        [~,mwu(:,2)] = getBHYqVls(mwu(:,1)',qv{2});
        
    otherwise
        % This uses the mafdr function obtained via Matlab. It seems to be
        % rather a little too generous, as all q values are lower than the
        % p values. Thus not sure it is correct...
        [~,ttp(:,2)] = mafdr(ttp(:,1));
        [~,mwu(:,2)] = mafdr(mwu(:,1));
end

% Ditch the empty values from the list...
bl = ~cellfun(@isempty,fcInfo(:,1));
fcInfo = fcInfo(bl,:);
mwu = mwu(bl,:);
ttp = ttp(bl,:);

% Save individual values for the mean/median thresholding
mw.mean   = ttp;
mw.median = mwu;

% Now export the data to a spreadsheet to include the pq values and the
% fold changes for the assigned lipids.
newName = [sP(1:end-4) '-FoldChangeData-' hID{1} '-' hID{2} '.txt'];
fID = fopen(newName,'w');
fprintf(fID,['Class\tLength\tDesaturations\tm/z\t',...
    'log2 mean FC\tlog2 median FC\t',...
    'ANOVA p\tANOVA q\t',...
    'MWU p\tMWU q\t']);
fprintf(fID,[hID{1} '-Mean\t' hID{2} '-Mean\t' hID{1} '-Median\t' hID{2} '-Median\t'...
    hID{1} '-STD\t' hID{2} '-STD\t' hID{1} '-25th\t' hID{1} '-75th\t'...
    hID{2} '-25th\t' hID{2} '-75th\n']); 
for n = 1:size(fcInfo,1)
    ii = fcInfo{n,6};
    
    % Determine the mean and median values for the two groups
    lipmean = [mean(gui.spec(indX,ii)) mean(gui.spec(indY,ii))];
    lipmed  = [median(gui.spec(indX,ii)) median(gui.spec(indY,ii))];
    
    % Standard deviations?
    lipstd = [std(gui.spec(indX,ii)) std(gui.spec(indY,ii))];
    lip25  = [prctile(gui.spec(indX,ii),25) prctile(gui.spec(indY,ii),25)];
    lip75  = [prctile(gui.spec(indX,ii),75) prctile(gui.spec(indY,ii),75)];

    % What is the m/z value of this variable?
    mzVal = sprintf('%0.4f',gui.op.cmz(ii));
    
    if isinf(fcInfo{n,4})
        fc1 = '';
    else
        fc1 = sprintf('%0.4f',fcInfo{n,4});
    end
    if isinf(fcInfo{n,5})
        fc2 = '';
    else
        fc2 = sprintf('%0.4f',fcInfo{n,5});
    end
    
    if length(fcInfo{n,7}) == 1
        wrClss = [fcInfo{n,1} '-' fcInfo{n,7}];
    else
        wrClss = fcInfo{n,1};
    end
        
    fprintf(fID,'%s\t%d\t%d\t%s\t%s\t%s\t%0.4E\t%0.4E\t%0.4E\t%0.4E\t',...
        wrClss,fcInfo{n,2},fcInfo{n,3},mzVal,fc1,fc2,...%fcInfo{n,4},fcInfo{n,5},...
        mw.mean(n,1),mw.mean(n,2),...
        mw.median(n,1),mw.median(n,2));
    
    % Here pritn the mean/median/std/iqr for each class
    fprintf(fID,'%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\n',...
        lipmean(1),lipmean(2),...
        lipmed(1),lipmed(2),...
        lipstd(1),lipstd(2),...
        lip25(1),lip75(1),...
        lip25(2),lip75(2));
    
end
fclose(fID);
disp('**** Fold change data saved to: ****');
disp(newName);
disp('**** ************************** ****');


% Set infinite fold change values to 0
infs = cellfun(@isinf,fcInfo(:,4:5));
for n = 1:size(infs,1)
    
    % Set the mean fold infinite changes to zero, and corresponding mwu p/q
    % values
    if infs(n,1) == 1
        fcInfo(n,4) = {0};
        mw.mean(n,:) = 1;
    end
    
    % Set the median fold infinite changes to zero, and the same for the
    % MWU-derived p/q values
    if infs(n,2) == 1
        fcInfo(n,5) = {0};
        %fcInfo(n,6) = {1};
        mw.median(n,:) = 1;
    end
end

% Ditch the insignificant entries if their q-value is too high, i.e.
% above the threshold...
keep = mw.median(:,2) <= qv{3};
fcInfoMEDIAN = fcInfo(keep,:);

keep = mw.mean(:,2) <= qv{3};
fcInfoMEAN = fcInfo(keep,:);

% Define the red/white/blue colormap
tr = linspace(0,1,50);
cols = ones(101,3);
cols(1:50,2) = tr;
cols(1:50,3) = tr;
cols(1:50,1)   = 1;
cols(52:101,1) = fliplr(tr);
cols(52:101,2) = fliplr(tr);
cols(52:101,3)   = 1;
cols = flipud(cols);

% Need to tell the figure which way the fold change has been calculated
txt = [hID{1} ' v ' hID{2}];

% DO MEANS FIRST
drawFoldChange(fcInfoMEAN  ,cols,'mean',txt,qv,true);
drawFoldChange(fcInfoMEAN  ,cols,'mean',txt,qv,false);
drawFoldChange(fcInfo,cols,'mean',txt,qv,false);

% NOW DO THE MEDIANS
drawFoldChange(fcInfoMEDIAN,cols,'median',txt,qv,true);
drawFoldChange(fcInfoMEDIAN,cols,'median',txt,qv,false);
drawFoldChange(fcInfo,cols,'median',txt,qv,false);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawFoldChange(fcInfo,cols,type,txt,pqt,labelFlag)
% Draw the fold changes nicely in one mega figure

if size(fcInfo,1) == 0
    return;
end

% Decide whether we do mean / median FCs
switch type
    case 'mean'
        ii = 4;
    case 'median'
        ii = 5;
end

% Create a figure and define three axes
fig.fig = figure('Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5],...
    'Color','white',...
    'Name',txt,...
    'Tag','meanmedianfoldchange'); 
%hold on;
fig.ax1 = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.05 0.05 0.25 0.9]);
fig.ax2 = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.35 0.05 0.25 0.9]);
fig.ax3 = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.65 0.05 0.25 0.9]);

uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0 0 0.2 0.02],...
    'Style','text',...    
    'HorizontalAlignment','center',...
    'String',['FDR = ' sprintf('%f',pqt{2}) ', q < ' sprintf('%f',pqt{3})],...
    'BackgroundColor','white');


drawIndividualHeatMap(fig.ax1,1,ii,fcInfo,cols,[],labelFlag);
drawIndividualHeatMap(fig.ax2,2,ii,fcInfo,cols,[],labelFlag);
drawIndividualHeatMap(fig.ax3,3,ii,fcInfo,cols,[],labelFlag);



% Finally add a single colour bar to the last of the images...
cb = colorbar;
set(cb,'Position',[0.925 0.05 0.025 0.9]);
set(cb,'FontSize',18);
ylabel(cb,['log_2 ' type ' fold change'],'FontSize',18,'FontWeight','bold');



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawIndividualHeatMap(ax,sortIndex,ii,fcInfo,cols,lims,labelFlag)
% PLace the individual heat map...

% If the limits are set, then remove values...
if ~isempty(lims)
    idx1 = cellfun(@lt,fcInfo(:,sortIndex),repmat({min(lims)},size(fcInfo,1),1));
    fcInfo(idx1,:) = [];
    
    idx2 = cellfun(@gt,fcInfo(:,sortIndex),repmat({max(lims)},size(fcInfo,1),1));
    fcInfo(idx2,:) = [];
end

% First we draw the classes in their place...
axes(ax);
hms = sortrows(fcInfo,[sortIndex ii]);
dat = cell2mat(hms(:,ii));

% Define color limits
clim = max(abs(dat(:)));

% Determine where to put the tick labels and lines across the classes
[lab,labInd,newInds] = tickLabDet(hms,size(dat,1),sortIndex);

% Draw the image with limits
imagesc(dat(:,1),[-clim clim]);

% Add dividing lines
for n = 1:numel(lab)
    line([0.5 1.5],[labInd(n)-0.5 labInd(n)-0.5],'LineWidth',1,'Color','black');
end
line([0.5 1.5],[size(dat,1)+0.5,size(dat,1)+0.5],'LineWidth',1,'Color','black');

% Boundary lines vertically
xlx = xlim;
line([xlx(1) xlx(1)],[0.5 size(dat,1)+0.5],'LineWidth',1,'Color','black');
line([xlx(2) xlx(2)],[0.5 size(dat,1)+0.5],'LineWidth',1,'Color','black');

% Format labels and ticks and stuff
set(gca,'XTick',[]);
set(gca,'YTick',newInds,'YTickLabel',lab,'TickLength',[0 0],'YDir','reverse');
colormap(cols);

set(gca,'FontSize',16,'FontWeight','bold');
ylim([0.499 size(dat,1)+0.499]);

% Now add the text labels for each of the lipids...
if labelFlag
    for n = 1:size(fcInfo,1)
        
        if length(hms{n,7}) == 1
            txtLab = [hms{n,1} '(' hms{n,7} '-' int2str(hms{n,2}) ':' int2str(hms{n,3}) ')'];
        else
            txtLab = [hms{n,1} '(' int2str(hms{n,2}) ':' int2str(hms{n,3}) ')'];
        end
        
        text(1,n,txtLab,...
            'HorizontalAlignment','center',...
            'FontSize',14);
    end
end



box off;

% Add in a title
switch sortIndex
    case 1
        axTit = 'Class';
    case 2
        axTit = 'Acyl chain length';
    case 3
        axTit = 'Desaturations';
end
title(axTit,'FontSize',18,'FontWeight','bold');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lab,labInd,newInds] = tickLabDet(hms,szDat,sortIndex)
% Determine tick label placement values

try
    [lab,labInd,~] = unique(hms(:,sortIndex));
catch
    [lab,labInd,~] = unique(cellfun(@unique,hms(:,sortIndex)));
end
newInds = zeros(numel(lab),1);
for n = 1:numel(lab)
    
    if n ~= numel(lab)
        if labInd(n+1)-labInd(n) == 1
            newInds(n,1) = labInd(n);
        else
            newInds(n,1) = mean(labInd(n:n+1));
        end
    else
        if szDat == labInd(n)
            newInds(n,1) = labInd(n);
        else
            newInds(n,1) = mean([labInd(n) szDat]);
        end
    end    
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g1,g2,lng,qv] = foldChangeDialog(hIDs)
% Ask the user which groups to use in the fold change

drpfig = figure('Units','normalized',...
    'Position',[0.45 0.40 0.1 0.2],...
    'Toolbar','none',...
    'Menubar','none',...
    'Name','Fold Change',...
    'Number','off');

h1 = uicontrol('Parent',drpfig,...
    'Style','popupmenu',...
    'Units','normalized',...
    'Position',[0.1 0.85 0.8 0.1],...
    'String',hIDs,...
    'Value',1);
    
h2 = uicontrol('Parent',drpfig,...
    'Style','popupmenu',...
    'Units','normalized',...
    'Position',[0.1 0.75 0.8 0.1],...
    'String',hIDs,...
    'Value',2);

h3 = uicontrol('Parent',drpfig,...
    'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.1 0.025 0.8 0.1],...
    'String','OK',...
    'Callback','uiresume(gcbf)');

q0 = uicontrol('Parent',drpfig,...
    'Style','popupmenu',...
    'Units','normalized',...
    'Position',[0.1 0.55 0.8 0.1],...
    'String',{'BHY'; 'MAFDR'},...
    'Value',1,...
    'Enable','off');

q1 = uicontrol('Parent',drpfig,...
    'Style','popupmenu',...
    'Units','normalized',...
    'Position',[0.1 0.45 0.8 0.1],...
    'String',{'FDR = 0.001'; 'FDR = 0.01'; 'FDR = 0.05'; 'FDR = 0.1'},...
    'Value',1);

q2 = uicontrol('Parent',drpfig,...
    'Style','popupmenu',...
    'Units','normalized',...
    'Position',[0.1 0.35 0.8 0.1],...
    'String',{'q < 0.001'; 'q < 0.01'; 'q < 0.05'; 'q < 0.1'},...
    'Value',1);


% Boxes for the determination of limits for length only
b1 = uicontrol('Parent',drpfig,...
    'Style','edit',...
    'Units','normalized',...
    'Position',[0.15 0.2 0.3 0.1],...
    'String','20');
b2 = uicontrol('Parent',drpfig,...
    'Style','edit',...
    'Units','normalized',...
    'Position',[0.55 0.2 0.3 0.1],...
    'String','60');


uiwait;

% Determine outputs
g1 = get(h1,'Value');
g2 = get(h2,'Value');
lng = [str2num(get(b1,'String')) str2num(get(b2,'String'))];

% Get the method...
switch get(q0,'Value')
    case 1
        qv{1} = 'bhy';
    case 2
        qv{1} = 'mafdr';
end

% Get the FDR rate for the BHY method
switch get(q1,'Value')
    case 1
        qv{2} = 0.001;
    case 2
        qv{2} = 0.01;
    case 3
        qv{2} = 0.05;
    case 4
        qv{2} = 0.1;
end

% Get the value for the threshold
switch get(q2,'Value')
    case 1
        qv{3} = 0.001;
    case 2
        qv{3} = 0.01;
    case 3
        qv{3} = 0.05;
    case 4
        qv{3} = 0.1;
end
    

close(drpfig);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
