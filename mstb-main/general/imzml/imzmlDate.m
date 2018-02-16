function [ tc ] = imzmlDate(file)
%imzmlDate - read in the date/time stamp

% Initial function
try
    
    tree = xmlread(file);

    iteration = 0;
    [~,~,res] = parseChildNodes(tree,iteration,[]);

    % Now we need to convert the datestr into a useable value...
    tc = datenum([res(1:10) res(12:end)],'yyyy-mm-ddHH:MM:SS');

catch
    tc = now;
    disp('Could not get date');
end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [children,iteration,res] = parseChildNodes(theNode,iteration,res)
% Recurse over node children

list = {'#text';'cvParam';'cvList';'sourceFile'};
children = [];
if theNode.hasChildNodes
    childNodes = theNode.getChildNodes;
    numChildNodes = childNodes.getLength;
    allocCell = cell(1,numChildNodes);
    children = struct('Name',allocCell,'Attributes',allocCell,...
                                 'Data',allocCell,'Children',allocCell);
    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        nname = char(theChild.getNodeName);
                
        % Things that we can ignore...
        if any(strcmp(list,nname))
            continue;
        end
        
        %disp([int2str(iteration) ' - ' nname]);
        
        % What we are really looking for...
        if strcmp(nname,'run')
            
            att = parseAttributes(theChild);
            
            % Check that we have the datestr
            fx = strcmp({att.Name},'startTimeStamp');
            res = att(fx).Value;
            return;
        end
           
        % Continue through to next child as necessary...
        if isempty(res)        
            % Just loop through each one if we don't find the one we want
            [~,iteration,res] = parseChildNodes(theChild,iteration+1,res);
        end
        
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function attributes = parseAttributes(theNode)
% Create attributes struct
attributes = [];
if theNode.hasAttributes
    theAttributes = theNode.getAttributes;
    numAttributes = theAttributes.getLength;
    allocCell = cell(1,numAttributes);
    attributes = struct('Name',allocCell,'Value',allocCell);
    for count = 1:numAttributes
        attrib = theAttributes.item(count-1);
        attributes(count).Name = char(attrib.getName);
        attributes(count).Value = char(attrib.getValue);
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nodeStruct = makeStructFromNode(theNode)

nodeStruct = struct('Name',char(theNode.getNodeName),...
    'Attributes',parseAttributes(theNode),'Data','',...
    'Children',parseChildNodes(theNode));

if any(strcmp(methods(theNode),'getData'))
   nodeStruct.Data = char(theNode.getData); 
else
    nodeStruct.Data = '';
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
