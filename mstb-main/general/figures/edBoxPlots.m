function [ output_args ] = edBoxPlots(mz,data,groups,lipid)
%luisaBoxPlots

% These are the colours of the groups
%cols = [19 179 232; 64 168 35; 149 23 191; 232 217 51]/256;

%cols = [227 23 168; 108 245 149] / 256; %cut
cols = [12 133 34; 237 29 19] / 256; % cutcoagcomb

% This is the order
order = [1 2];

% Get the pertinent information
[mzList,fileSave,mult] = getValues(lipid);

% Find the m/z values
idx = mzFind(mz,mzList,2);
find(idx)

% Trim the list of variables...
tmp = data(:,idx);

%tmp(tmp(:,1) > 5e5,1) = NaN;

% Run the boxplot function
jsmBoxPlotMulti(tmp,groups,round(mz(idx)*100)/100,...
    'Colours',cols,...
    'Order',order,...
    'Orientation','vertical',...
    'Legend',true,...
    'Mult',mult,...
    'Labels',{'Normal','Cancer'});%,'Cancer'});


% Save the graph...
graphFormat(fileSave,'png')


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mzList,fileSave,mult] = getValues(lipid)

% These are the mz values to be found
switch lower(lipid)
    
    case 'pa'
        mzList = [699.50; 701.52; 723.50; 725.52; 747.50];
        mult = [5 1 5 5 5];
        
    case 'ps'
        mzList = [760.52 796.55 836.54 842.60 882.53]';
        mult = [2 1 2 4 5];
        
    case 'pe'
        mzList = [720.50 762.51 796.59 816.56 818.58]';
        mult = [2 1 3 5 5];
        
    case 'pi'
        mzList = [819.53 847.57 865.57 889.57]';
        mult = [1 4 2 1];
        
    case 'cer'
        mzList = [654.56 656.58 682.59]';
        mult = [6 6 1];
        
    case 'esj-cancer'
        mzList = [687.50 716.52 744.55 766.54 772.58]';
        mult = [2 2 1 2 2];
        mult = [1 1 1 1 1];
        
    case 'esj-normal'
        mzList = [865.70 891.72 893.73 917.73 921.76];
        mult = [3 2 1 3 3];
        %mult = [1 1 1 1 1];
end

fileSave = ['/Users/jmckenzi/Desktop/Boxplots-' upper(lipid)];

end