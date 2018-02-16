function isoGen
%isoGen - generate m/z values and an accurate isotopic distribution for a
%given molecule
%
%
% James McKenzie, 2016

% Draw the window
[fig] = drawWindow;

set(fig.go,'Callback',{@calcFunc,fig});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = drawWindow
% Draw the window

fig.fig = figure('Units','normalized',...
    'Position',[0.25 0.25 0.5 0.4]);

% We'll have an axes on the right hand side
fig.ax = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.55 0.05 0.4 0.9]);

% And on the left we'll have boxes and options for elements, and other
% options
fig.form = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.005 0.9 0.2 0.095],...
    'Style','edit',...
    'String','C6H6',...
    'FontSize',20,...
    'BackgroundColor',[0.9 0.9 0.9]);

% Button to calculate the composition / distrubution / etc
fig.go = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.215 0.9 0.2 0.095],...
    'Style','pushbutton',...
    'String','Calculate',...
    'FontSize',20,...
    'BackgroundColor',[0.9 0.9 0.9]);



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calcFunc(~,~,fig)
% Function that does the calculations and stuff

% Get the formula
form = get(fig.form,'String');

% Determine what it is
elem = parseFormula(upper(form));

% Calculate distribution
[dist] = isotopicCalculator(elem.qty);

% Plot it!
stem(1:numel(dist),dist,...
    'LineWidth',10,...
    'MarkerSize',1);


% Calculate the masses of the elements
totMass = elem.mono .* elem.qty;

% Empty table
tabDat = cell(size(elem.mono,2)+1,3);
for n = 1:size(elem.mono,2)
    tabDat(n,1) = elem.symb(n);
    if elem.qty(n) > 0
        tabDat(n,2) = {elem.qty(n)};
        tabDat(n,3) = {totMass(n)};
    end
end
tabDat(n+1,1) = {'Total'};
tabDat(n+1,3) = {sum(totMass)};

% Trim out empty entries
fx = cellfun(@isempty,tabDat(:,3));
tabDat = tabDat(~fx,:);

xlim([0.5 numel(dist)+0.5]);
ylim([0 105]);

tabDat

totMass'

% Total mass in loooong format
disp(['Total mass = ' sprintf('%0.6f',tabDat{end,end})]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%