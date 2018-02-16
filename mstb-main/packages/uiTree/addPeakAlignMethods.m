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
treeModel.insertNodeInto(parentNode,node,node.getChildCount());
maxMZShiftNode = uitreenode('v0','param', 'max mz shift: 1', [], 0);
%jmenu.add('max mz shift');
treeModel.insertNodeInto(maxMZShiftNode,RSPANode,RSPANode.getChildCount());
minPeakWidth   = uitreenode('v0','param', 'min peak width: 1', [], 0);
treeModel.insertNodeInto(minPeakWidth,RSPANode,RSPANode.getChildCount());
recursion      = uitreenode('v0','param', 'recursion: 1', [], 0);
treeModel.insertNodeInto(recursion,RSPANode,RSPANode.getChildCount());
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