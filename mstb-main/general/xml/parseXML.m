function [solv,flow,mz] = parseXML(filename)
% PARSEXML Convert XML file to a MATLAB structure.
try
    tree = xmlread(filename);
catch
    error('Failed to read XML file %s.',filename);
end

% Recurse over child nodes. This could run into problems
% with very deeply nested trees.
%try
    [theStruct,quit] = parseChildNodes(tree,false);
%catch
%    error('Unable to parse XML file %s.',filename);
%end


% Find the bits of info that I actually want...
s1 = theStruct(1).Children(14).Children(2).Children(6).Children(2).Children(8).Attributes(3).Value;
solv = theStruct(1).Children(14).Children(2).Children(6).Children(2).Children(8).Attributes(4).Value;

s2 = theStruct(1).Children(14).Children(2).Children(6).Children(2).Children(6).Attributes(3).Value;
flow = theStruct(1).Children(14).Children(2).Children(6).Children(2).Children(6).Attributes(7).Value;

s3 = theStruct(1).Children(6).Children(6).Children(6).Attributes(4).Value;
br1 = strfind(s3,'[');
br2 = strfind(s3,']');
mz = s3(br1+1:br2-1);

%disp([s1 ' - ' solv]);
%disp([s2 ' - ' flow]);
%disp(['mz' ' - ' mz]);

end


% ----- Local function PARSECHILDNODES -----
function [children,quit] = parseChildNodes(theNode,quit)
% Recurse over node children.

if quit
    children = [];
    return
end

children = [];
if theNode.hasChildNodes
    childNodes = theNode.getChildNodes;
    numChildNodes = childNodes.getLength;
    allocCell = cell(1, numChildNodes);
    
    children = struct(             ...
        'Name', allocCell, 'Attributes', allocCell,    ...
        'Data', allocCell, 'Children', allocCell);
    
    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        try
            [children(count),quit] = makeStructFromNode(theChild,quit);
        catch
            quit = true;
            return
        end
        
        %if strcmp(children(count).Attributes,'STOP')
        %    return
        %end
    end
end

if quit
    return
end

end

% ----- Local function MAKESTRUCTFROMNODE -----
function [nodeStruct,quit] = makeStructFromNode(theNode,quit)
% Create structure of node info.

if quit
    nodeStruct = [];
    return
end

[atts,quit1] = parseAttributes(theNode,quit);
[chil,quit2] = parseChildNodes(theNode,quit);
if quit1 || quit2
    quit = true;
else
    quit = false;
end

nodeStruct = struct(                        ...
    'Name', char(theNode.getNodeName),       ...
    'Attributes', atts,  ...
    'Data', '',                              ...
    'Children', chil);

% Decide if it is here where we stop...
if strcmp(nodeStruct.Attributes,'STOP')
    return
end

if any(strcmp(methods(theNode), 'getData'))
    nodeStruct.Data = char(theNode.getData);
else
    nodeStruct.Data = '';
end

end

% ----- Local function PARSEATTRIBUTES -----
function [attributes,quit] = parseAttributes(theNode,quit)
% Create attributes structure.

if quit
    return
end

attributes = [];
if theNode.hasAttributes
    theAttributes = theNode.getAttributes;
    numAttributes = theAttributes.getLength;
    allocCell = cell(1, numAttributes);
    attributes = struct('Name', allocCell, 'Value', ...
        allocCell);
    
    for count = 1:numAttributes
        attrib = theAttributes.item(count-1);
        attributes(count).Name = char(attrib.getName);
        attributes(count).Value = char(attrib.getValue);
    end

    %    meta = cell(1,4);
    %    if size(attributes,2) >= 3
    %        if strcmp(attributes(3).Value,'solvent')
    %            meta{1,1} = attributes(4).Value;
    %
    %        elseif strcmp(attributes(3).Value,'filter string')
    %            meta{1,2} = attributes(4).Value;
    %
    %        elseif strcmp(attributes(3).Value,'solvent flowrate')
    %            meta{1,3} = attributes(4).Value;
    %
    %        elseif strcmp(attributes(3).Value,'spray voltage')
    %            meta{1,4} = attributes(7).Value;
    %
    %        end
    %        %meta
    %    end
    
    
end

% Here is where we would consider stopping...
try
    if strcmp(attributes(1).Value,'MS:1000544')
        %disp('Time to stop');
        quit = true;
        return
    end
catch err    
end

end