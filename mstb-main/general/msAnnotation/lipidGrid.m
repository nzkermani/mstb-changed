function [ tabs ] = lipidGrid(op,ass,fold)
%lipidGrid - uses the op and assignments to create a range of lipid maps
%for length .v. desaturations etc...

warning on all;

f0 = findobj('Name','Lipid Grid');
close(f0);

op.XPeaks = full(op.XPeaks);
op.XPeaksLog = full(op.XPeaksLog);

% Remove all NaN/Inf values from the data...
nnn = isnan(op.XPeaks);
op.XPeaks(nnn) = 0;
iii = isinf(op.XPeaks);
op.XPeaks(iii) = 0;
sum(nnn(:))
sum(iii(:))

% Make the storage elements
[tabs] = makeTables(op,ass);

% Add in the information
[tabs] = popTables(op,ass,tabs);

% Then you draw a window for displaying the tables and performing stocsy...
[fig] = drawWindow(op,tabs);

% Run the initial 'show all' populate...
%subplotPopulate([],[],'median',fig,tabs)

% Enable callbacks
set(fig.drawMeanFC,  'Callback',{@drawSubplots,fig,tabs,op,'mean'});
set(fig.drawMedianFC,'Callback',{@drawSubplots,fig,tabs,op,'median'});
set(fig.drawSTOCSY,  'Callback',{@drawSubplots,fig,tabs,op,'STOCSY'});
set(fig.drawboxp,    'Callback',{@drawBoxPlot,op});
set(fig.selMZ,       'Callback',{@selectSTOCSY,op});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tabs] = makeTables(op,ass)
% Make and fill the tables for info storage

% Define the parameters for searching
lipLen = 30:44;
fasLen =  10:30;
limSat =  0:10;

% So need to do this for 'every' lipid class (PA,PE,PG,PI,PS etc), so need
% to find all the unique subclasses (ignore classes!)
lipCls = unique(ass.db.c2);
numLip = numel(lipCls);
tabs = struct('class',[],'length',[],'desat',[],...
    'varNum',[],...
    'mean',[],...
    'median',[],...
    'stoCor',[],...
    'stoCov',[]);

% Determine some more key variables
%numE     = numel(op.cmz);
histIDs  = unique(op.histID);
numHist  = numel(histIDs);

% Make empty tables for the structure - there will be a lot of information
% to be stored and gathererd...
for n = 1:numLip
    
    % Class name
    tabs(n).class = lipCls{n};
    
    % Some subclasses will have specific lengths limits...
    switch lipCls{n}
        case 'FA'
            tabs(n).length = fasLen;
        otherwise
            tabs(n).length = lipLen;
    end
    
    % Currently all the same desat count
    tabs(n).desat = limSat;   
    
    % Short lenghts
    tmpLen = numel(tabs(n).length);
    tmpSat = numel(tabs(n).desat);
    
    % Empty table for mean/median intensities
    tabs(n).mean    = zeros(tmpLen,tmpSat,numHist);
    tabs(n).median  = zeros(size(tabs(n).mean));
    
    % Somewhere for variable numbers...
    tabs(n).varNum = zeros(tmpLen,tmpSat);
    
    % STOCSY - for storage perhaps...
    %...eventually...    
       
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tabs] = popTables(op,ass,tabs)
% Here we fill the tables with stuff

% Define necessary parameters
numE    = numel(op.cmz);
lipCls  = unique(ass.db.c2);
histIDs = unique(op.histID);
numHist = numel(histIDs);

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
        
        warning('What are you doing about P- and O- lipids?');
        if ~isempty(cho) && ~isempty(cho) && ~isempty(chl) %%&& ~skip
            
            % How many carbons?
            nCarb  = tmp(cho(1)+1:chl(1)-1);

            % This might include a P- or an O-...            
            if strcmp(nCarb(1:2),'P-') || strcmp(nCarb(1:2),'O-')
                nCarb = nCarb(3:end);
            end
                        
            % Should be straight forward
            nDesat = tmp(chl(1)+1:chc(1)-1);
            
            % Convert to numbers
            nCarb = str2double(nCarb);
            nDesat= str2double(nDesat);
        
            % To which class does this belong?
            annoCls = ass.annoCls{n,1};

            % What is the true mz of the assignment?
            %fx  = find(ass.list(:,1) == n);
            %ffx = ass.list(fx,2);         %#ok<FNDSB>
            %trueMZ = ass.db.mz(ffx);

            % ppm difference
            %ppmDev = 1e6 * ( gui.data(1).mz(1,n) - trueMZ) / trueMZ;
            %ppmDev = sprintf('%0.2f',ppmDev);
            %trueMZ = sprintf('%0.4f',trueMZ);
            
            % So now that we know everything, we need to put it in the
            % correct place of the table structure.
            chc = strfind(annoCls,',');
            annoCls2 = annoCls(chc(1)+1:end);
            fz = strcmp(lipCls,annoCls2);
            fz = find(fz == 1);

            % First the row (length)
            [fx] = find(tabs(fz).length == nCarb);

            % Now the sat (column)
            [fy] = find(tabs(fz).desat == nDesat);
            
            % Add the variable number to the table
            tabs(fz).varNum(fx,fy) = n;

            % Take the STOCSY info from the structure
            %tabs(fz).stoCov(fx,fy) = gui.stocsy(end).res.CovXY(1,n);
            %tabs(fz).stoCor(fx,fy) = gui.stocsy(end).res.CCXY( 1,n);
            
            % Now place various stuff in the various tables...
            for r = 1:numHist            
                
                % Which are the correct entries?
                [tx] = strcmp(op.histID,histIDs{r});

                % Calculate the mean/median of this tissue class
                tabs(fz).median(fx,fy,r) = median(op.XPeaks(tx,n));
                tabs(fz).  mean(fx,fy,r) =   mean(op.XPeaks(tx,n));

            end
        else
            % Couldn't add it in...            
        end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = drawWindow(op,tabs)
% Draw a window for displaying the results...

% Get the names of the lipid classes
lipClass = cell(size(tabs,2),1);
for n = 1:size(tabs,2)
    lipClass{n,1} = tabs(n).class;
end

% Draw the figure
fig.fig = figure('Name','Lipid Grid',...
    'Number','off',...
    'Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5],...
    'Toolbar','figure');

% Need to add in some boxes to specify the histIDs for calculation of the
% lipid fold changes and stuff...
histIDs = unique(op.histID);
%numH = numel(histIDs);

fig.fc1 = uicontrol('Parent',fig.fig,...
    'Style','popupmenu',...
    'Units','pixels',...
    'Position',[0 0 100 30],...
    'String',histIDs,...
    'Value',1);

fig.fc2 = uicontrol('Parent',fig.fig,...
    'Style','popupmenu',...
    'Units','pixels',...
    'Position',[100 0 100 30],...
    'String',histIDs,...
    'Value',2);

% Button for STOCSY!
fig.selMZ = uicontrol('Parent',fig.fig,...
    'Style','pushbutton',...
    'Units','pixels',...
    'Position',[5 35 90 30],...
    'String','STOCSY',...
    'UserData',1);

% How about being able to select which lipid classes to include?
fig.inc = uicontrol('Parent',fig.fig,...
    'Style','listbox',...
    'Units','normalized',...
    'Position',[0.0 0.8 0.035 0.2],...
    'String',lipClass,...
    'Value',1:numel(lipClass),...%[5 7 9 10 11 12],...
    'Max',size(lipClass,1));

% Button for drawing...
fig.drawMeanFC = uicontrol('Parent',fig.fig,...
    'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.0 0.75 0.035 0.05],...
    'String','MeaFC');

fig.drawMedianFC = uicontrol('Parent',fig.fig,...
    'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.0 0.70 0.035 0.05],...
    'String','MedFC');

fig.drawSTOCSY = uicontrol('Parent',fig.fig,...
    'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.0 0.65 0.035 0.05],...
    'String','STOCSY');

fig.drawboxp = uicontrol('Parent',fig.fig,...
    'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.0 0.60 0.035 0.05],...
    'String','Box!');

% Tick boxes to display the variable numbers / variable masses, which will
% be useful for seeing what has been assigned, and whether it is 'correct'.
fig.boxNum = uicontrol('Parent',fig.fig,...
    'Style','checkbox',...
    'Units','pixels',...
    'Position',[0 100 100 30],...
    'String','Var #',...
    'Value',0);

fig.boxMZ = uicontrol('Parent',fig.fig,...
    'Style','checkbox',...
    'Units','pixels',...
    'Position',[0 130 100 30],...
    'String','Var m/z',...
    'Value',1);


% Now draw the subplots and stuff! Mean fold change by default
drawSubplots([],[],fig,tabs,op,'mean');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawSubplots(src,event,fig,tabs,op,method)
% Add in a couple of subplots...

% Find out which ones are selected!
vals = get(fig.inc,'Value');

% And their names?
strn = get(fig.inc,'String');
%strn = strn(vals);

% How many figures to draw?
numFig = numel(vals);

numC = ceil(sqrt(numFig));
numR = ceil(numFig / numC);

% Delete all existing axes (in the current figure!
f0 = findall(gcf,'Type','axes');
delete(f0);

% Start with the subplots...
sp = zeros(numFig+1,1);
for n = 1:numFig
    sp(n,1) = subplot(numR+1,numC,n+numC);    
end

% Axes at the top to display the spectral variables...
sp(numFig+1,1) = subplot(numR+1,numC,1:numC);

% Need to have a generic colorbar
scb = colorbar;
set(scb,'Position',[0.925 0.1 0.025 0.8]);
caxis([-1 +1]);

% Set as user data
gui.axHand  = sp;
gui.scb     = scb;

set(fig.inc,'UserData',gui);

% Now actually add in the stuff...
subplotPopulate(src,event,method,fig,tabs,op)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function subplotPopulate(src,event,method,fig,tabs,op)
% Populate the subplots with something useful

cmap = jet(100);
cmap(52,:) = [1 1 1];
cmap(49,:) = [1 1 1];
cmap(50,:) = [1 1 1];
cmap(51,:) = [1 1 1];

% Need to get the decision of what to plot!
vals = get(fig.inc,'Value');

% And their names?
strn = get(fig.inc,'String');
%strn = strn(vals);

% How many figures to draw?
numFig = numel(vals);
        
% Place to stores imgs
img = struct('name',[],'img',[],'min',[],'max',[]);

% What to plot?
switch method    
    case {'mean','median'}
        % Much of the method is essentially the same, so we'll do that
        % instead of duplicating code twice...
        
        % Need to get the histIDs to be fold changed
        grpA = get(fig.fc1, 'Value');
        grpB = get(fig.fc2, 'Value');

        for n = 1:numel(vals)
            
            % Need to say which one this corresponds to...
            i = vals(n);
            
            % Here we just prepare the image
            switch method
                case 'mean'
                    tmp = log2(tabs(i).mean(:,:,grpA) ./ ...
                        tabs(i).mean(:,:,grpB));
                case 'median'
                    tmp = log2(tabs(i).median(:,:,grpA) ./ ...
                        tabs(i).median(:,:,grpB));
            end
            
            % Set nans infs to zero
            tmp(isnan(tmp)) = 0;
            tmp(isinf(tmp)) = 0;
            
            % Add in
            img(n).img = tmp;
            
            img(n).min = min(tmp(:));
            img(n).max = max(tmp(:));
        end
        
        % Now find the axes limits...
        minI = min(vertcat(img.min));
        maxI = max(vertcat(img.max));
        absI = max([abs(minI) abs(maxI)]);
        clim = [-absI +absI];
            
        yLab = ['log_2 ' method ' fold change'];
        
    case 'STOCSY'
        % Calcualte only the values that are necessary, thus it will be as
        % quick as possible...
        
        % To which variable are we correlating?
        stoInd = get(fig.selMZ,'UserData');
        if isempty(stoInd)
            warning('Select an m/z value first!');
            return
        end
        
        % Loop through the selected things
        for n = 1:numel(vals)
            
            % Need to say to which this corresponds
            i = vals(n);
            
            % Get the indices of the variables for this group...
            inds = tabs(i).varNum;
            szI  = size(inds);
            crs  = zeros(szI(1),szI(2));
            
            for j = 1:szI(1)
                for k = 1:szI(2)
                    
                    if inds(j,k) ~= 0
                        % Then correlate...
                        %crs(j,k) = j*k;
                        [rho,pval] = corr(op.XPeaks(:,stoInd),op.XPeaks(:,inds(j,k)));
                        if pval <= 0.05
                            crs(j,k) = rho;
                        end
                    else
                        % Variable not found, so no correlation possible
                    end
                end
            end
            
            % Set NaNs / INFs to zero
            crs(isnan(crs)) = 0;
            crs(isinf(crs)) = 0;
            
            % Add in to the structure
            img(n).img = crs;
            img(n).min = -1;
            img(n).max = +1;            
        end
        
        % Now find axes limits
        absI = 1;
        clim = [-1 +1];
        
        yLab = ['STOCSY: ' get(fig.selMZ,'String')];
        
    otherwise
        % There is no otherwise
end

% Are we going to include the text strings of variable names / mzs? Need to
% get the values of the tick boxes to see...
varNum = get(fig.boxNum,'Value');
varMZ  = get(fig.boxMZ, 'Value');

% Now that the images have been prepared and are ready to go, let's plot
% them...
tmp = get(fig.inc,'UserData');
for n = 1:numel(tmp.axHand)-1
    
    % Set the current axes
    axes(tmp.axHand(n,1)); %#ok<LAXES>
    cla reset
    hold on;
    
    % Finally plot...
    i = vals(n);
    if clim(1) == 0 && clim(2) == 0
        imagesc(tabs(i).desat,tabs(i).length,img(n).img);
    else
        imagesc(tabs(i).desat,tabs(i).length,img(n).img,clim);
    end
    colormap(cmap);    
    
    % Add in the variable annotations / numbers
    if varNum == 1 || varMZ == 1
        for k = 1:numel(tabs(i).length)
            for j = 1:numel(tabs(i).desat)
                                
                % Here we prepare the string...and add it to the axes                
                if tabs(i).varNum(k,j) == 0
                    tmpTxt = '.';
                    
                    text(tabs(i).desat(j),tabs(i).length(k),...
                        tmpTxt,...
                        'Units','data',...
                        'Clipping','on',...
                        'HorizontalAlignment','Center',...
                        'Color','black');
                    
                else
                    if varNum == 1 && varMZ == 1
                        tmpTxt = [int2str(tabs(i).varNum(k,j)) ', '...
                            sprintf('%0.4f',op.cmz(tabs(i).varNum(k,j)))];
                    elseif varNum == 1
                        tmpTxt = int2str(tabs(i).varNum(k,j));
                    elseif varMZ == 1
                        tmpTxt = sprintf('%0.4f',op.cmz(tabs(i).varNum(k,j)));
                    end
                    
                    % Color choice? Determine the colour choice according
                    % to its colour: if green/blue in the middle it is
                    % better black, if dark red/blue better white...
                    pct = 100 * abs(img(n).img(k,j)) / max(clim);

                    if pct > 50
                        tCol = 'white';
                    else
                        tCol = 'black';
                    end
                    

%                     if img(n).img(k,j) == 0
%                         tCol = 'black';
%                     else                        
%                         tCol = 'white';
%                     end
                                
                    text(tabs(i).desat(j),tabs(i).length(k),...
                        tmpTxt,...
                        'Units','data',...
                        'Clipping','on',...
                        'HorizontalAlignment','Center',...
                        'Color',tCol,...
                        'FontWeight','bold');

                        
                    
                    
                    
                end
            end
        end
    end
               
    % The addition of the text means that we have to specify the axes
    % limits to keep the images tight to the axes
    xlim([min(tabs(i).desat)-0.5 max(tabs(i).desat)+0.5]);
    ylim([min(tabs(i).length)-0.5 max(tabs(i).length)+0.5]);

    set(gca,'YDir','reverse');
    
    % Add the lipid class title at the top...
    title(tabs(i).class, 'FontSize',14,'FontWeight','bold');
end

% Prepare the cbar tick labels
yticks = [-absI:absI/5:0 absI/5:absI/5:absI];
yticks = round(yticks * 100) / 100;

% Add the label etc...
ylabel(tmp.scb,...
    yLab,...
    'FontSize',14,...
    'FontWeight','bold');

caxis(clim);
colormap(cmap);
set(tmp.scb,'YTickLabel',yticks)


% Now need to plot the data... but i tihnk just the mean of the histIDs
% rather than everything else, as there are too many data points.
unqHist = unique(op.histID);

axes(tmp.axHand(end));
hold on;

% Define some colours...
cols = jet(numel(unqHist));

unqFile = unique(op.pID);
hh = zeros(numel(unqHist),1);
legLabs = cell(numel(unqHist),1);

% Create a 'dummy' mz vector with inserted zeros so that the ploting is a
% little better...
zeroSpacing = 0.01;


for n = 1:numel(unqFile)
    [fx] = strcmp(op.pID,unqFile{n});
    
    
    % Now for each of the histIDs in this file...
    fileHist = unique(op.histID(fx));
    for r = 1:numel(fileHist)        
        [fy] = strcmp(op.histID,fileHist{r});        
        fz = (fx + fy) == 2;
        
        % Determine which colour it should be...
        ci = strcmp(unqHist,fileHist{r});
        
        switch method
            case 'mean'
                [xp,yp] = insertZeros(op.cmz,mean(op.XPeaks(fz,:),1),zeroSpacing);                
            case 'median'
                [xp,yp] = insertZeros(op.cmz,mean(op.XPeaks(fz,:),1),zeroSpacing);                
            otherwise
                [xp,yp] = insertZeros(op.cmz,mean(op.XPeaks(fz,:),1),zeroSpacing);                
        end
        hh(ci) = plot(xp,yp,'Color',cols(ci,:));
        legLabs{ci,1} = fileHist{r};
    end
    
    
end
legend(hh,legLabs);
xlim([min(xp) max(xp)]);

    



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectSTOCSY(src,event,op)
% Run stocsy on a single m/z value

% First start by gathering a single m/z value using ginput...
[mz,~] = ginput(1);

if mz < 50
    warning('You''ve clicked in the wrong axes most probably.');
    return
end

% Convert the selected mz to the closest variable
px = abs(op.cmz - mz);
[~,px] = min(px);

% Just update the STOCSY mz value here. Will trigger another function to
% run the actual analysis.
set(src,...
    'String',['m/z = ' sprintf('%0.4f',op.cmz(px))],...
    'UserData',px);

% Now run STOCSY-like analysis using this selected peak... But feasibly
% don't need to do it for all variables? Just the ones that feature in the
% lipid analysis? Likely to make it quicker...

% Start by running through the data to get the variable numbers, then
% correlate and shove the 



%tmp = full(op.XPeaks);
%crr = corr(tmp,tmp(:,px));

%figure; imagesc(crr);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawBoxPlot(src,event,op)
% Draw a box plot for a single variable

% First start by gathering a single m/z value using ginput...
[mz,~] = ginput(1);

if mz < 50
    warning('You''ve clicked in the wrong axes most probably.');
    return
end

% Convert the selected mz to the closest variable
px = abs(op.cmz - mz);
[~,px] = min(px);

hfig = figure('Name','Box Plot');
hold on;

% Determine x positions according to histID
[unqN,~,unqI] = unique(op.histID);

% Add a little wiggle to the data...
wigval = 0.5;
wiggle = rand(numel(unqI),1) * wigval;
wiggle = wiggle - wigval/2;
xval = unqI + wiggle;

% Use this data...
ydata = op.XPeaksLog(:,px);

% Scatter the points
scatter(xval,ydata,50,'blue','o');

% Add the mean and median values for each group...
for n = 1:numel(unqN)
    
    % Indices
    fx = find(unqI == n);
    
    % Mean
    mn = mean(ydata(fx,:));
    
    % Median
    md = median(ydata(fx,:));
    
    % Draw on plot...
    hnd(1) = line([n-wigval/2 n+wigval/2],[mn mn],...
        'LineWidth',4,'Color','red');
    
    hnd(2) = line([n-wigval/2 n+wigval/2],[md md],...
        'LineWidth',4,'Color','magenta');
end

xlim([min(unqI)-0.5 max(unqI)+0.5]);
legend(hnd,{'Mean','Median'},'Location','NorthEastOutside',...
    'Orientation','vertical');

title(['m/z = ' sprintf('%0.4f',op.cmz(px))],'FontSize',18);
ylabel('Logged intensity','FontSize',14);

set(gca,'XTick',1:numel(unqN),'XTickLabel',unqN)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%