function uiSetupPreprocParams(preproc,parentIndex,nodeIndex,jtree,tree,treeModel,rootNode,node,figHand)
% uiSetupPreprocParams - customization of uitree pre-processing parameters 
%  Input:   preproc     - tree like structure of all parameters
%           parentIndex - parent node index
%           nodeIndex   - node index
% Author:   Kirill A. Veselkov, Imperial College 2014.
% Improved: James S. McKenzie, Imperial College 2014.

methodid  = preproc{parentIndex}.methodnames{nodeIndex};
defparams = preproc{parentIndex}.defparams{nodeIndex};
if isfield(preproc{parentIndex},'options') ...
        && (length(preproc{parentIndex}.options)>=nodeIndex)
    options = preproc{parentIndex}.options{nodeIndex};
else
    options = [];
end

hFigParam = figure('MenuBar','none',...
    'ToolBar','none',...
    'units','pixels',...
    'Visible','off'); 
updateFigTitleAndIconMS(hFigParam, methodid,'MSINavigatorLogo.png'); 
drawnow; 

nParams         = length(defparams)/2;
columnname      = cell(1,nParams); 
columnformat    = cell(1,nParams); 
columneditable  = ones(1,nParams)==1;
for index = 1:2:nParams*2
    iParam             = (index+1)/2;
    columnname{iParam} = defparams{index};
    if isempty(options) || length(options)<iParam
        columnformat{iParam}   = 'char';
    else
        if ~isempty(options{iParam})
            columnformat{iParam} = options{iParam};
        end
    end
end

% Make the table allowing you to change the setttings
hParamTable = uitable('Units','pixels',...
    'Position',[0 0 1 1],...
    'Data', defparams(2:2:end),...
    'ColumnName',columnname,...
    'ColumnFormat', columnformat,...
    'ColumnEditable', columneditable,...
    'RowName',[],...
    'Tag','uitreenodedata');%,...
    %'CellEditCallback',{@defParamEditCallBack,preproc,jtree,'none'});
    
% Adjust positions
tablePstns  = get(hParamTable,'Extent');   % Required uitable size
if nParams<4
    tablePstns(3) = tablePstns(3);%./nParams; % minimum number of parameters
end
%set(hParamTable,'Position',tablePstns); 

% Setup "OK" and "Make as defaults" push buttons
PB1Pstns    = tablePstns;
PB1Pstns(3) = tablePstns(3);
PB1Pstns(4) = (tablePstns(4));
PB1Pstns(2) = 0;

tablePstns(2) = tablePstns(2)+PB1Pstns(4);

uicontrol(hFigParam,...
    'Style','PushButton',...
    'Units','pixels',...
    'String','Set',...
    'Position',PB1Pstns,...
    'Tag','Set',...
    'Callback',{@defParamEditCallBack,preproc,jtree,tree,hParamTable,treeModel,rootNode,node,figHand});

%PB2Pstns    = PB1Pstns;
%PB2Pstns(1) = PB1Pstns(1)+PB1Pstns(3);

% uicontrol(hFigParam,...
%     'Style','PushButton',...
%     'Units','pixels',...
%     'String','Make as defaults',...
%     'Tag','Make as defaults',...
%     'Position',PB2Pstns,...
%     'Callback',{@defParamEditCallBack,preproc,jtree,'set'});


set(hParamTable,'Position',tablePstns); 
Pstns = get(hFigParam,'Position');

Pstns(4) = (PB1Pstns(4) + tablePstns(4)); 
Pstns(3) = tablePstns(3); 
set(hFigParam,'Position', Pstns,'Visible','on');
uitreeIndexPath = [parentIndex,nodeIndex];
guidata(hFigParam,uitreeIndexPath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function defParamEditCallBack(hFig,ignore,preproc,jtree,tree,tabHand,treeModel,rootNode,node,figHand)
% Callback function - get the values from the cell(s) in the table, and
% 'write' them to the options structure (which i might have to make...)

% Get table data...
tabDat = get(tabHand,'Data');

% These are the node indices...
i = guidata(hFig);

% Update in this loop...
for n = 1:size(tabDat,2)    
    j = n*2;
    preproc{i(1)}.defparams{i(2)}(j) = tabDat(n);    
end

targetNode = node.getParent.getParent;
targetNode.removeAllChildren(); 
treeModel.reload(targetNode);
% update tree
addnodetopreprocroot( preproc{i(1)},treeModel,...
                                rootNode, targetNode);
jtree.treeDidChange();
tree.expand(targetNode);
close(get(tabHand,'Parent'));
                            
% Perhaps consider changing the defparams permanently here - not sure
% exactly how. Perhaps get Kirill to finish it.

% Need to update the tree handle into the user data thing of the msDB,
% which hopefully will solve the problem of not being able to read in the
% 'new' parameters (currently only reading in the original ones)
%f0 = findobj('Tag','uiTreePre');
%set(f0,'UserData',tree);        

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%