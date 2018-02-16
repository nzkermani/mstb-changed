function [ output_args ] = annaBoxPlots(mz,data,groups,type)
%luisaBoxPlots

% These are the colours of the groups
%cols = [19 179 232; 64 168 35; 149 23 191; 232 217 51; 100 100 100]/256;
cols = [46 58 230; 250 143 2]/256;

% This is the order
order = [];

% Get the pertinent information
[mzList,fileSave,mult,histID] = getValues(type);

% Find the m/z values
idx = mzFind(mz,mzList,2);

% Histological information
numG = numel(histID);
include = zeros(size(data,1),1);
for n = 1:numG
    tmp = strcmp(groups,histID{n});
    include = include + tmp;
end
inc = include > 0;
    

% Trim the list of variables...
tmp = data(inc,idx);

% Run the boxplot function
jsmBoxPlotMulti(tmp,groups(inc,:),round(mz(idx)*100)/100,...
    'Colours',cols,...
    'Order',order,...
    'Orientation','vertical',...
    'Legend',true,...
    'Mult',mult);

% Save the graph...
graphFormat(fileSave,'png')


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mzList,fileSave,mult,histID] = getValues(group)

% These are the mz values to be found
switch lower(group)
    
    case 'ald-nash'
        mzList = [682.5962;688.4945;770.5726;792.5345;909.5526];
        mult = [1 1 1 1 1];
        histID = {'ALD Nodules';'NASH Nodules'};
    
    case 'aih-pbc'
        mzList = [684.6095;817.5055;819.5216;865.5055];
        mult = [6 3 1 5];
        histID = {'AIH Nodules';'PBC Nodules'};
        
    case 'fn'
        mzList = [750.5465;788.5845;861.5525;887.5725;909.5526];
        mult = [5 5 1 3 5];
        histID = {'Fibrosis';'Nodules'};

end

fileSave = ['/Users/jmckenzi/Desktop/Anna-BP-' upper(group)];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
