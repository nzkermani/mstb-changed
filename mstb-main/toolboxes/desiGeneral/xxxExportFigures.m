function xxxExportFigures(src,event,fig)
% xxxExportFigures - from the toolbox, we copy the old figures into new
% ones.

% Full method taken from: https://uk.mathworks.com/matlabcentral/
% answers/93226-how-do-i-copy-the-axes-from-an-existing-gui-
% into-a-new-figure-so-that-the-new-figure-does-not-have-a
titleSize = 18;

% Get the guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% desi/desiPosNeg, need to determine appropriate axes
[axs,txt,defs] = determineAxes(fig,dpn.mode);

% Ask the user which one we want to export
[sel,~] = listdlg('PromptString','Select axes',...
    'SelectionMode','single',...
    'ListString',txt,...
    'ListSize',[160 160]);
if isempty(sel)
    return
end

% Create new figure and copy to it
figNew = figure('Units','normalized',...
    'Position',defs{sel}{1});
copyobj(axs(sel),figNew);

% Get the axes handle
g = findobj('Parent',figNew,'Type','axes');

% Set the new position
switch txt{sel}
    case {'Optical','Positive','Negative','Fused'}
        set(g,'Position',[0.05 0.05 0.9 0.9]);
    otherwise
        set(g,'Position',[0.1 0.1 0.8 0.8]);
end

% Set axes size
set(g,'FontSize',14,'FontWeight','bold');

% Change the title (if present)
getTitleText(g,titleSize);
%title('');

% Perhaps we can delete the patches? (manual un/comment)
ptch = findobj('Parent',gca,'Type','patch');
delete(ptch);

% Depending on the type of image, change the various properties...
switch defs{sel}{2}
    
    case 'opt'
        
    case 'ms'
        colormap(redbluecmap);
        axis square
        
    case 'load'
        
    case 'scatter'
        
    otherwise
        disp('No otherwise');
end

% Increase the title's font size

% Add the m/z label
%xlabel('m/z  ','FontSize',14,'FontWeight','bold','FontAngle','italic');


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [axs,txt,defs] = determineAxes(fig,mode)
% Determine the axes to show in the plot

% Determine view mode
flag = strcmp(get(fig.tb.layout,'State'),'off');

% Define a few image sizes and stuff in here...
szs = {[0.25 0.25 0.2 0.35],[0.25 0.25 0.2 0.25]};


% Get the axes for each one
if flag && strcmp(mode,'single')
    % Simple view / single ion mode
    axs = [fig.ax.opt(1) fig.ax.ms1(1)];
    txt = {'Optical','Mass Spec'};
    defs{1} = {szs{1},'opt'};
    defs{2} = {szs{1},'ms'};

elseif flag && strcmp(mode,'dual')
    % Simple view / dual ion mode
    axs = [fig.ax.opt(1) fig.ax.ms1(1) fig.ax.ms2(1) fig.ax.fu(1)];
    txt = {'Optical','Positive','Negative','Fused'};
    defs{1} = {szs{1},'opt'};
    defs{2} = {szs{1},'ms'};
    defs{3} = {szs{1},'ms'};
    defs{4} = {szs{1},'ms'};

elseif ~flag && strcmp(mode,'single')
    % Complex view / single ion mode
    axs = [fig.ax.opt(1) fig.ax.ms1(1) fig.ax.sp(1) fig.ax.mv(1)];
    txt = {'Optical','Mass Spec','Loadings','Scores'};
    defs{1} = {szs{1},'opt'};
    defs{2} = {szs{1},'ms'};
    defs{3} = {szs{1},'load'};
    defs{4} = {szs{1},'scatter'};

elseif ~flag && strcmp(mode,'dual')
    % Complex view / dual ion mode
    axs = [fig.ax.opt(1) fig.ax.ms1(1) fig.ax.ms2(1) fig.ax.fu(1) fig.ax.sp fig.ax.mv fig.ax.mv];
    txt = {'Optical','Positive','Negative','Fused','Top Right','Bottom Right','BR-Small'};
    defs{1} = {szs{1},'opt'};
    defs{2} = {szs{1},'ms'};
    defs{3} = {szs{1},'ms'};
    defs{4} = {szs{1},'ms'};
    defs{5} = {szs{2},'load'};
    defs{6} = {szs{1},'scatter'};
    defs{7} = {szs{2},'scatter'};

else
    error('Nothing good here');
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function getTitleText(g,fS)

%title(g,'');
%g = findobj(gca,'Type','patch'); delete(g);

%return

try
    titTxt = get(g.Title,'String');
    
    if strcmp(titTxt(1:3),'m/z')
        com1 = strfind(titTxt,',');
        titTxt = titTxt(1:com1(1)-1);
    end
    
    % See if it is a Latex string
    ltx = strfind(titTxt,'\fontsize{');
    if ~isempty(ltx)
        braces = strfind(titTxt,'}');        
        idx = [ltx(1)+length('\fontsize{') braces(1)-1];        
        titTxt = [titTxt(1:idx(1)-1) int2str(fS) titTxt(braces(1):end)];        
    end
        
    
    title(g,titTxt,'FontSize',fS,'FontWeight','bold');

catch
    disp('Cannot set title');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%