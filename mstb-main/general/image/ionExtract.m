function ionExtract(varargin)
%ionExtract - master function that extracts ions, via a graphical user
%interface...

% Input arguments and stuff...
[opt] = readArgsData(varargin);

% First we draw a window that will contain all the interesting bits.  We'll
% need a series of controls down one side, and an image axes down the rest
% of it
[fig] = drawWindow([opt.fP opt.fN]);

% Set the callback functions here...
set(fig.mode,'Callback',{@changeMode});
set(fig.file,'Callback',{@fileSelect});
set(fig.go,  'Callback',{@doExtract});
set(fig.spectrum,'Callback',{@plotAvgSpec});
set(fig.figExport,'Callback',{@figureExport});

% Save the guidata
ie.fig = fig;
ie.opt = opt;
guidata(ie.fig.fig,ie);

% Run initial functions...
changeMode(ie.fig.fig,[]);
set(ie.fig.mass,'String',num2str(ie.opt.ions));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = readArgsData(argsin)
% Read the arguments and then the data if it wasn't passed

% These are the defaults
opts.fP     = '/Users/jmckenzi/DB/ColorectalDummy/A30/';
opts.fN     = 'A30centreS475um.imzML';
opts.inst   = 'orbi';

% Need to define the adduct lists of the various possibilities...
adds = struct('mode',[],'list',[],'name',[],'elem',[]);
adds(1).mode = 'neg';
adds(2).mode = 'pos';
[adds(1).list,adds(1).name,adds(1).elem] = adductLists('n');
[adds(2).list,adds(2).name,adds(2).elem] = adductLists('p');
opts.adds = adds;

% Pre-define some ions here
%opts.ions = [885.55; 886.55];
opts.ions = [280.2038;290.1518;296.0467;296.3079;298.178;300.1937;...
    310.2872;316.2038;318.2195;320.2285;320.2351;320.2351;322.2508;...
    324.1937;324.2664;325.2627;326.1649;326.2093;326.2722;328.1522;...
    328.1805;328.225;328.2845;330.1679;332.1736;332.1988;332.1988;...
    334.2144;334.2144;336.23;336.2301;336.2301;336.2301;336.2301];

% There are too many here for testing
% opts.ions = [280.2038;290.1518;296.0467;296.3079;298.178;300.1937;...
%     310.2872;316.2038;318.2195;320.2285;320.2351;320.2351;322.2508;...
%     324.1937;324.2664;325.2627;326.1649;326.2093;326.2722;328.1522;...
%     328.1805;328.225;328.2845;330.1679;332.1736;332.1988;332.1988;...
%     334.2144;334.2144;336.23;336.2301;336.2301;336.2301;336.2301;...
%     338.2457;338.2457;338.2457;338.2457;340.1886;340.2548;340.2614;...
%     342.2042;342.2406;342.277;344.2351;348.1937;348.2301;350.2093;...
%     350.2093;350.2093;350.2093;350.2093;350.2093;350.2457;352.1805;...
%     352.225;352.225;352.225;352.225;352.225;352.225;352.225;352.225;...
%     354.2406;354.2406;354.277;356.2497;356.2563;358.2653;362.2457;...
%     363.2773;364.1886;364.2614;365.3086;366.2042;366.2042;366.2042;...
%     366.2042;366.2406;367.3086;368.2199;368.2199;368.2199;368.2199;...
%     368.2563;370.1756;370.2355;370.2355;370.2355;372.1937;372.1937;...
%     372.2446;372.2512;372.2512;374.2602;374.2602;378.277;379.2723;...
%     380.1754;380.2563;381.2879;382.1992;382.2719;382.2719;384.2148;...
%     386.2093;386.2305;387.1352;387.241;388.225;390.2018;390.2406;...
%     393.2879;395.2672;396.1704;396.2512;396.2876;396.2876;397.2828;...
%     397.2828;401.3075;402.2254;408.2512;410.2668;413.1508;414.2254;...
%     415.1665;418.3083;424.1249;425.2777;426.2618;427.2934;428.2774;...
%     430.1328;430.2719;438.1809;439.2392;440.1966;441.1821;441.2549;...
%     444.2148;444.2723;446.2305;453.2185;454.2567;455.2342;455.2342;...
%     457.2498;468.2723;469.2134;472.107;481.2498;484.1461;494.2451;...
%     496.1864;496.2607;498.202;500.2177;504.2359;512.1813;512.2556;...
%     514.197;514.2713;532.1281;538.197;540.2126;556.2075;558.2232;...
%     568.2818;598.2181;600.2337;625.3033;641.2982;643.3139];

% Ditch any replicated values.
opts.ions = unique(opts.ions);

return

% Run through each pair
nArgs = length(argsin);
for i = 1:2:nArgs
    if strcmpi('mz',argsin{i})
        opts.mz = argsin{i+1};
        
    elseif strcmpi('folder',argsin{i})
        opts.fP = argsin{i+1};

    elseif strcmpi('file',argsin{i})
        opts.fN = argsin{i+1};
    
    end
end

if ~exist([opts.fP opts.fN],'file')
    error('File does not exist');
end

if numel(opts.mz) ~= 2
    error('Give two mz values over a range');
end

if opts.mz(2) < opts.mz(1)
    error('Specify mz range properly');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = drawWindow(file)
% Draw the window and add all the bits to it

f0 = findobj('Tag','ionExtract');
close(f0);

fig.fig = figure('Name','ionExtract',...
    'Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5],...
    'Tag','ionExtract');


% Main axes
fig.ax = uipanel('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.4 0 0.6 1],...
    'BackgroundColor',[0 0 0]);

% Now a panel in which to place all the controls
fig.pan = uipanel('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0 0.3 0.4 0.7],...
    'BackgroundColor',[0.8 0.8 0.8]);

% Panel bits and bobs
fig.mass = uicontrol('Parent',fig.pan,...
    'Units','normalized',...
    'Position',[0 0 0.24 0.9],...
    'Style','edit',...
    'Min',1,...
    'Max',10,...
    'FontSize',16);

fig.mode = uicontrol('Parent',fig.pan,...
    'Units','normalized',...
    'Position',[0.25 0.8 0.24 0.1],...
    'Style','popupmenu',...
    'String',{'Negative','Positive'},...
    'Value',1);

fig.ions = uicontrol('Parent',fig.pan,...
    'Units','normalized',...
    'Position',[0.25 0.25 0.24 0.55],...
    'Style','listbox',...
    'FontSize',14,...
    'Min',1,...
    'Max',10);

uicontrol('Parent',fig.pan,...
    'Units','normalized',...
    'Position',[0.25 0.07 0.24 0.055],...
    'Style','text',...
    'String','+/- ppm',...
    'FontSize',16,...
    'BackgroundColor',[0.8 0.8 0.8]);

fig.ppm = uicontrol('Parent',fig.pan,...
    'Units','normalized',...
    'Position',[0.25 0 0.24 0.05],...
    'Style','edit',...
    'String','10',...
    'FontSize',16);

if ~exist(file,'file')
    file = 'Select File';
end
fig.file = uicontrol('Parent',fig.pan,...
    'Units','normalized',...
    'Position',[0 0.95 1 0.05],...
    'Style','pushbutton',...
    'String',file,...
    'FontSize',14);

% A button to go...
fig.go = uicontrol('Parent',fig.pan,...
    'Units','normalized',...
    'Position',[0.75 0 0.25 0.1],...
    'Style','pushbutton',...
    'String','Go');

% Secondard panel for results section...
fig.pan2 = uipanel('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0 0 0.4 0.3],...
    'BackgroundColor',[0.8 0.8 0.8]);

% Text string to say which ion mode was analysed
fig.modeDone = uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0 0.9 0.24 0.1],...
    'Style','text',...
    'String','None',...
    'FontSize',16,...
    'BackgroundColor',[0.8 0.8 0.8]);

% Box for the masses to be displayed in
fig.plot = uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0 0 0.24 0.8],...
    'Style','listbox',...
    'FontSize',14,...
    'Min',1,...
    'Max',1);

% Global colour limit
fig.clim = uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.25 0.7 0.49 0.1],...
    'Style','checkbox',...
    'String','Intensity Limit?',...
    'FontSize',14,...
    'BackgroundColor',[0.8 0.8 0.8]);

% Log transform box
fig.log = uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.25 0.6 0.49 0.1],...
    'Style','checkbox',...
    'String','Log Intensities?',...
    'FontSize',14,...
    'BackgroundColor',[0.8 0.8 0.8]);

% Log transform box
fig.waters = uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.25 0.5 0.49 0.1],...
    'Style','checkbox',...
    'String','Waters Line Removal?',...
    'FontSize',14,...
    'BackgroundColor',[0.8 0.8 0.8]);

% Smoothing option
fig.smooth = uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.25 0.4 0.14 0.075],...
    'Style','popupmenu',...
    'String',{'1','2','3','4','5','10'},...
    'Value',1);

uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.4 0.375 0.14 0.1],...
    'Style','text',...
    'String','Smoothing',...
    'FontSize',14,...
    'BackgroundColor',[0.8 0.8 0.8]);

% Button for the average/sum spectrum
fig.spectrum = uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.75 0.6 0.24 0.2],...
    'Style','pushbutton',...
    'String','Spectrum');

fig.figExport = uicontrol('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.75 0.4 0.24 0.2],...
    'Style','pushbutton',...
    'String','Figure Export');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeMode(src,event)
% Change ion mode from neg / pos

% Guidata
ie = guidata(src);

% Selection box
val = get(ie.fig.mode,'Value');
txt = get(ie.fig.mode,'String');

% Get the text labels of the adducts
switch txt{val}    
    case 'Negative'
        adds = ie.opt.adds(1).name;        
    case 'Positive'
        adds = ie.opt.adds(2).name;
    otherwise
        % There is no otherwise
end

% Update the box of adducts
set(ie.fig.ions,'String',adds,'Value',1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fileSelect(src,event)

[fN,fP,~] = uigetfile({'*.imzML';'*.DAT';'*.dat'});

if isnumeric(fN)
    set(src,'String','Select File');
    return
end

set(src,'String',[fP fN]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doExtract(src,event)

% Guidata
ie = guidata(src);

% Determine the various bits of information that we need in order to
% proceed
file = get(ie.fig.file,'String');
if ~exist(file,'file')
    disp('No file');
    return
end

% Ion mode
mode = get(ie.fig.mode,'String');
tmp  = get(ie.fig.mode,'Value');
mode = mode{tmp};

% Selected adducts... (the indices will do for now...
adds = get(ie.fig.ions,'Value');

% Now the mass values to be analysed
mass = get(ie.fig.mass,'String');
mass = str2num(mass); %#ok<ST2NM>
numMass = numel(mass);

ppmTol = get(ie.fig.ppm,'String');
ppmTol = str2num(ppmTol); %#ok<ST2NM>


% Now we need to do the extraction for these masses --> ions in the file.
% This requires us to convert each mass into the series of ion m/z values,
% and these are extracted from the file in turn. It would be excellent to
% display the different ion forms alongside each other, so save the results
% into a structure, with the images stored as an [m x n x p] datacube
allMZ = cell(numMass,4);
for n = 1:numMass
    
    % Convert the mass to the selected ion forms
    switch lower(mode)        
        case 'negative'
            val = 1;           
        case 'positive'
            val = 2;
    end
    
    % Calculate values of m/z
    [mzgrid] = mass2mz(mass(n),ie.opt.adds(val));
        
    % Trim out the ones we don't want
    mzVals = mzgrid(adds)';
    mzName = ie.opt.adds(val).name(adds);
    
    % Save to a generic structure
    allMZ{n,1} = mass(n);
    allMZ{n,2} = mzName;
    allMZ{n,3} = mzVals;

end

% Determine what kind of file (imzML or Waters) that we are trying to 
% extract from. This will only work for Waters RAW files, not those by any
% other manufacturer...
extnDot = strfind(file,'.');
extn = lower(file(extnDot(end)+1:end))

% Now with a single list of m/z values, we open the imzML file and start to
% do the business
ionsTmp = vertcat(allMZ{:,3});
switch extn
    case 'imzml'
        [extr,spectrum] = imzmlIonExtract(file,ionsTmp,ppmTol);
        
    case 'dat'
        % Going to be a little temperamental I expect...
        try
            extnSlash = strfind(file,filesep);
            file = file(1:extnSlash(end)-1);
            [extr,spectrum] = extractWatersRaw(file,ionsTmp,ppmTol);
        catch err
            err
            error('Possibly the wrong kind of file / platform / day...');
        end
end

% Now we need to consider separating the ions in order to be able to
% display these images properly
i = 0;
for n = 1:numel(adds):numel(ionsTmp)
    i = i + 1;
    st = n;
    fn = n + numel(adds) - 1;    
    allMZ{i,4} = extr(:,:,st:fn);
end
    
% Save the results to the gui structure
ie.extr = allMZ;
ie.spec = spectrum;

% Change the various things for plotting...
set(ie.fig.plot,...
    'String',num2str(mass),...
    'Value',1,...
    'Callback',{@doPlots});

% Tell us what mode we've extracted
set(ie.fig.modeDone,'String',mode);

% Update the guidata...
guidata(src,ie);


% ... then call the plotting function
doPlots(ie.fig.fig,[]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doPlots(src,event)

% Get guidata
ie = guidata(src);

figure(ie.fig.fig);

% Delete existing axes?
f0 = findobj('Tag','iesp2');
delete(f0);

% Determine which mass should be plotted
mass = get(ie.fig.plot,'String');
tmp  = get(ie.fig.plot,'Value');

% Determine the number of ions that have been extracted
numI = size(ie.extr{tmp,4},3);

% What is the best arrangement?
optRC = sqrt(numI);
if isinteger(optRC)
    arr = [optRC optRC];
else
    tmp2 = ceil(optRC);
    optRC = ceil(numI / tmp2);
    arr = [optRC tmp2];
end

% This is the ion data
name = ie.extr{tmp,2};
list = ie.extr{tmp,3};
ions = ie.extr{tmp,4};

% We need to decide if we should get rid of the lines at the top of the
% images, as these are plagueing Waters data. This should be done before
% all efforts at colour scaling / normalisation / transformation
watersLine = get(ie.fig.waters,'Value')
if watersLine
    ions = ions(6:end,:,:);
end

% Need to see if we should log the data?
doLog = get(ie.fig.log,'Value');
if doLog
    ions = log(ions + min(ions(ions > 0)));
    logMessage = 'Based on log transformed intensities';
else
    logMessage = 'Based on raw intensities';
end

% What about smoothing?
smTmp = get(ie.fig.smooth,'Value');
xxTmp = get(ie.fig.smooth,'String');
smVal = str2double(xxTmp(smTmp));
if smVal == 1
    % Do nothing, i.e. no smoothing required
    smoMessage = 'No smoothing applied';
else
    % Need to perform smoothing on each of the images. First define the
    % filter to be used
    fh = fspecial('average',smVal);
    
    % Now loop through each of the images
    for n = 1:numI                
        ions(:,:,n) = filter2(fh,ions(:,:,n));        
    end
    smoMessage = ['Smoothing to level ' int2str(smVal)];
end


% ppm tolerances
ppmTol = get(ie.fig.ppm,'String');
ppmTol = str2num(ppmTol); %#ok<ST2NM>

% Need to decide on a colour map scale limit
colLims = [min(ions(:)) max(ions(:))];
cLim = get(ie.fig.clim,'Value');

% Can we calculate the correlation coefficient betweeen all of these ion
% images in order to provide an idea of the relation between them?
reshaped = reshape(ions,[size(ions,1)*size(ions,2) size(ions,3)]);
[allCorr,pval] = corr(reshaped);
mask = pval > 0.05;
%mask(eye(size(allCorr,1)) == 1) = 1;
allCorr(mask) = 0;

% Set diagonal elements to 1
mask = eye(size(allCorr,1));
allCorr(mask == 1) = 1;

% Round the correlations a little
allCorr = round(allCorr * 100) / 100;

name2 = strrep(name,'-','_');
try
    tabCorr = array2table(allCorr,'RowNames',name,'VariableNames',name2);
catch
    tabCorr = array2table(allCorr,'RowNames',name);
end

% Display results in the matlab environment
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['Correlation Results for M = ' mass(tmp,:)]);
disp(logMessage);
disp(smoMessage);
disp(char(10));
disp(tabCorr);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(char(10));

% Draw each ion
for n = 1:numI
    
    % Subplot here
    spax(n,1) = subplot(arr(1),arr(2),n,...
        'Parent',ie.fig.ax);
    
    % This actually plots the ions...
    imagesc(ions(:,:,n));
    
    % Global colour limit application
    if cLim
        caxis(colLims);
    end
    
    % Need to have a title for each of the images in order to know what it
    % is showing us.
    ppmDiff = ppmTol * list(n) / 1e6;
    titTxt = ['m/z = ' sprintf('%0.3f',list(n)) ' ± ' ...
        sprintf('%0.3f',ppmDiff) ', ' name{n}];
    title(titTxt,'FontSize',14,'Color',[0.5 0.5 0.5]);
    
    
end

% Black to white via red colour map
colormap(hot);

% Not wholly necessary, but good to have nonetheless
linkaxes(spax,'xy');

% Axes formatting
set(spax,...
    'XTickLabel',[],...
    'YTickLabel',[],...
    'XTick',[],...
    'YTick',[],...
    'XColor',[0.2 0.2 0.2],...
    'YColor',[0.2 0.2 0.2],...
    'Tag','iesp2');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotAvgSpec(src,event)
% Plot the average spectrum, if available

ie = guidata(src);

if ~isfield(ie,'spec')
    disp('No spectrum');
    return
end

figure;
stem(ie.spec(:,1),ie.spec(:,2));
xlabel('m/z','FontSize',18,'FontWeight','bold');
ylabel('Total intensity','FontSize',18,'FontWeight','bold');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figureExport(src,event)
% Ask the user for a folder and then export all the relevant figures...
% Could take a while...

ie = guidata(src);

if ~isfield(ie,'extr')
    disp('No ions extracted');
    return
end

% Ask the user for a folder if James' default doesn't exist
tmpDir = 'C:\Users\jmckenzi\Desktop\IonExport\';
if ~exist(tmpDir,'dir')
    tmpDir = uigetdir;
    if length(tmpDir) < 2
        disp('No Directory Selected');
        return
    end
    if ~strcmp(tmpDir(end),filesep)
        tmpDir = [tmpDir filesep];
    end
end
disp(tmpDir);

wb = waitbar(0,'Exporting - don''t click elsewhere!','WindowStyle','modal');

% Here we draw the single figure that gets reused every time. The invert
% hard copy property allows us to save the black coloured background as is,
% rather than it defaulting to white
expFig = figure('Name','tmpfig',...
    'Units','pixels',...
    'Position',[1200 0 600 600],...
    'Color','black',...
    'Visible','off',...
    'InvertHardCopy','off');
expHand = axes('Parent',expFig);
%expCBar = colorbar;
colormap(hot);

% PPM tolerances
ppmTol = get(ie.fig.ppm,'String');
ppmTol = str2num(ppmTol); %#ok<ST2NM>


% Now we need to begin the mega automation of doing every single image...
numMass = size(ie.extr,1);
for n = 1:numMass
    
    % Decide on the folder name, which requires removing . and replacing
    % them with -
    fName = sprintf('%0.4f',ie.extr{n,1});
    fName = strrep(fName,'.','-');
    
    newFold = [tmpDir fName filesep];
    if ~exist(newFold,'dir')
        mkdir(newFold);
    end
    
    % How many ions in this image?
    numIons = size(ie.extr{n,2},1);
    for r = 1:numIons
        
        tmpDat = ie.extr{n,4}(:,:,r);
        tmpIon = ie.extr{n,3}(r);
        tmpNam = ie.extr{n,2}{r};
        
        % With this information make a figure...
        imagesc(tmpDat,'Parent',expHand);
        axis square
        box on;
        c = colorbar;
        
        set(expHand,...
            'XTickLabel',[],...
            'YTickLabel',[],...
            'XTick',[],...
            'YTick',[],...
            'XColor',[0.2 0.2 0.2],...
            'YColor',[0.2 0.2 0.2],...
            'Tag','ghj');
        
        set(c,...
            'XColor',[0.5 0.5 0.5],...
            'YColor',[0.5 0.5 0.5],...
            'FontSize',14);

        ppmDiff = ppmTol * tmpIon / 1e6;
        titTxt = ['m/z = ' sprintf('%0.3f',tmpIon) ' ± ' ...
            sprintf('%0.3f',ppmDiff) ', ' tmpNam];
        title(titTxt,'FontSize',14,'Color',[0.5 0.5 0.5]);
        
        % Generate the appropriate file name...
        newName = [fName '--' tmpNam '-' sprintf('%d',ppmTol) 'ppm'];
        
        % Save the figure...
        print(gcf, '-dpng', '-r200',[newFold newName]);
        
    end
    
    waitbar(n/numMass,wb,'Exporting - don''t click elsewhere!');   
    
    
end
    
close(expFig);

delete(wb);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
