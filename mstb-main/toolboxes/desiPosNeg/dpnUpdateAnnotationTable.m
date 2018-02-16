function dpnUpdateAnnotationTable(dpn,man)

try
    unqCols = cell2mat(dpn.anno(:,2));
    [~,b,~] = unique(unqCols);

    % Update the table...
    tabDat = dpn.anno(b,[4 5]);
    set(man.tab,'Data',tabDat,...
        'FontSize',16,...
        'ColumnWidth',{100 100});
catch
    % Presumably no data to be shown
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%