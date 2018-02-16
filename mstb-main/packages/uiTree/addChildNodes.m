function addChildNodes(treeModel,parentNode,childNodeName,selected,varargin)
% Profile Alignment 
% Put your algorithm here: Recursive segment wise peak alignment

if strcmp(selected,'selected')
    [I,map]   = checkedIcon;
    javaImage = im2java(I,map);
elseif strcmp(selected,'unselected')
    [I,map]   = uncheckedIcon;
    javaImage = im2java(I,map);
end
node       = uitreenode('v0',selected, childNodeName, [], 0);
node.setIcon(javaImage);
treeModel.insertNodeInto(node,parentNode,parentNode.getChildCount());
if ~isempty(varargin)
    nParams = length(varargin);
    for i = 1:2:nParams
        iParam = varargin(i:i+1);
        if isnumeric(iParam{2})
            iParam{2} = num2str(iParam{2});
        end
        iParamNode = uitreenode('v0',[childNodeName ': param'],...
            [iParam{1},': ',iParam{2}], [], 0);
        treeModel.insertNodeInto(iParamNode,node,...
            node.getChildCount());
    end
end
end

function [I,map] = checkedIcon()
I = uint8(...
    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0;
    2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,1;
    2,2,2,2,2,2,2,2,2,2,2,2,0,2,3,1;
    2,2,1,1,1,1,1,1,1,1,1,0,2,2,3,1;
    2,2,1,1,1,1,1,1,1,1,0,1,2,2,3,1;
    2,2,1,1,1,1,1,1,1,0,1,1,2,2,3,1;
    2,2,1,1,1,1,1,1,0,0,1,1,2,2,3,1;
    2,2,1,0,0,1,1,0,0,1,1,1,2,2,3,1;
    2,2,1,1,0,0,0,0,1,1,1,1,2,2,3,1;
    2,2,1,1,0,0,0,0,1,1,1,1,2,2,3,1;
    2,2,1,1,1,0,0,1,1,1,1,1,2,2,3,1;
    2,2,1,1,1,0,1,1,1,1,1,1,2,2,3,1;
    2,2,1,1,1,1,1,1,1,1,1,1,2,2,3,1;
    2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,1;
    2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,1;
    1,3,3,3,3,3,3,3,3,3,3,3,3,3,3,1]);
map = [0.023529,0.4902,0;
    1,1,1;
    0,0,0;
    0.50196,0.50196,0.50196;
    0.50196,0.50196,0.50196;
    0,0,0;
    0,0,0;
    0,0,0];
end

function [I,map] = uncheckedIcon()
I = uint8(...
    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1;
    2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,1;
    2,2,1,1,1,1,1,1,1,1,1,1,2,2,3,1;
    2,2,1,1,1,1,1,1,1,1,1,1,2,2,3,1;
    2,2,1,1,1,1,1,1,1,1,1,1,2,2,3,1;
    2,2,1,1,1,1,1,1,1,1,1,1,2,2,3,1;
    2,2,1,1,1,1,1,1,1,1,1,1,2,2,3,1;
    2,2,1,1,1,1,1,1,1,1,1,1,2,2,3,1;
    2,2,1,1,1,1,1,1,1,1,1,1,2,2,3,1;
    2,2,1,1,1,1,1,1,1,1,1,1,2,2,3,1;
    2,2,1,1,1,1,1,1,1,1,1,1,2,2,3,1;
    2,2,1,1,1,1,1,1,1,1,1,1,2,2,3,1;
    2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,1;
    2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,1;
    1,3,3,3,3,3,3,3,3,3,3,3,3,3,3,1]);
map = ...
    [0.023529,0.4902,0;
    1,1,1;
    0,0,0;
    0.50196,0.50196,0.50196;
    0.50196,0.50196,0.50196;
    0,0,0;
    0,0,0;
    0,0,0];
end