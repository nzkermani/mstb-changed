function [ pre ] = readUItree( src,event,pre,tree )
%readUItree - get the values of the stuff within it...
% James S. McKenzie, Imperial College London, 2014.

% Don't touch these at all!
tree1 = tree.tree;
root  = tree1.getPathForLocation(1,1);
node  = root.getLastPathComponent;

% Number of children nodes
numN = node.getChildCount;

% An empty cell for storing the stuff
%pre = cell(numN,1);

for n = 1:numN
    
    % This is the actual node
    ntmp = node.getChildAt(n-1);
        
    % Is it ticked?
    tick = ntmp.getValue;
    if strcmp(tick,'selected')
        pre{n}.do = 1;
    else
        pre{n}.do = 0;
    end
    
    % This is the name
    name = deHyper(ntmp.getName);
    pre{n}.name = name;
    %disp(['>>> ' char(name)]);

    % How many leaves does it have?
    child  = ntmp.getChildCount;
    % don't update pre-processing variable if the node is not selected
    if child == 0
        continue;
    end
    for r = 1:child
        
        % Name of child
        cname = deHyper(ntmp.getChildAt(r-1).getName);
        pre{n}.methodnames{r} = cname;
        %disp([char(9) '>>> ' char(cname)]);
        
        % Selected?
        ctick = ntmp.getChildAt(r-1).getValue;
        if strcmp(ctick,'selected')
            pre{n}.selected(1,r) = 1;
        else
            pre{n}.selected(1,r) = 0;
        end
        
        % Further subchildren-based nodes
        c2 = ntmp.getChildAt(r-1).getChildCount;  
        
        for s = 1:c2            
            cnode = ntmp.getChildAt(r-1).getChildAt(s-1);            
            c2 = deHyper(cnode.getName);
            [par,val] = unPair(c2);
            %disp([char(9) char(9) '>>> ' char(c2)]);
            if s == 1
                params = {par,val};
            else
                params = {params{:} par val};
            end
        end
        if c2 > 0
            pre{n}.defparams{r} = params;
            clear params;
        else
            pre{n}.defparams{r} = {};               
        end
    end
end


% This is the format into which to place the actual options
% preproc{i}.do               = 1;
% preproc{i}.name             = 'Source of Raw Data'; 
% preproc{i}.methodnames{1}   = 'Annotations of Histologist with same sample MMC prediction';
% preproc{i}.defparams{1}     = {'Probability Threshold / %',99};
% preproc{i}.methodnames{2}   = 'Annotations of Histologist';
% preproc{i}.selected         = [0, 1];

%assignin('base','pre',pre)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [deh] = deHyper(txt)
% De-hypertext a string assuming that the last '>' is the end of html

txt = char(txt);
grt = strfind(txt,'>');

if ~isempty(grt)    
    deh = txt(grt(numel(grt))+1:end);
else
    deh = txt;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [par,val] = unPair(txt)
% Take a string and split the stuff before and after a ': '

txt = char(txt);
sep = strfind(txt,':');

if ~isempty(sep)
    par = txt(1:sep(1)-1);
    val = txt(sep(1)+2:end);
else
    par = txt;
    val = [];
    disp('shome mistake shurely');
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
