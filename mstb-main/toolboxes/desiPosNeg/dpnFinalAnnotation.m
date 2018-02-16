function dpnFinalAnnotation(~,~,fig,man)
%dpnFinalAnnotation - close the window and update what needs to be updated

% Guidata
dpn = guidata(fig.fig);
if ~isfield(dpn,'anno')
    %delete(man.fig);
    return
end

% Sub-window table
tab = get(man.tab,'Data');
numT = size(tab,1);
tab2 = cell(numT,2);
for n = 1:numT
    
    colgen = tab{n,1};
    
    f1 = strfind(colgen,'<TD>');
    f2 = strfind(colgen,'</TD>');
    
    % This is the ID of the annotation...
    tab2{n,1} = str2num(colgen(f1+4:f2-1)); %#ok<ST2NM>
    
    % This is the text that has been provided for the annotation
    tab2{n,2} = tab{n,2};
    
end

% Annotations - set them to match what is in tab2
anno = dpn.anno;
numA = size(anno,1);
for n = 1:numA    
    usd = anno{n,2};    
    fx = cell2mat(tab2(:,1)) == usd;    
    anno{n,5} = tab2{fx,2};
end

% Update the guidata...
dpn.anno = anno;
guidata(fig.fig,dpn);

% In the end we need to ensure that we close the figure
%delete(man.fig);

end