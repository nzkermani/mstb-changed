function [ grp ] = statsObservationLabels(meta,str,val)
%statsObservationLabels - format the meta labels for each observation

% Single or multiple?
if numel(val) > 1
    
    % First group
    grp = meta.(str{val(1)});
    if strcmp(str{val(1)},'date')
        grp = statsDateConvert(grp);
    end

    if isnumeric(grp)
        grpX = cell(size(grp));
        for r = 1:size(grp,1)
            grpX{r,1} = int2str(grp(r));
        end
        grp = grpX;
    end
    
    % Subsequent groups
    for r = 2:numel(val)
        grp2 = meta.(str{val(r)});
        if strcmp(str{val(r)},'date')
            grp2 = statsDateConvert(grp2);
        end
        grp = cat(2,grp,grp2);
    end
    
    % Combine
    grp = classMany2One(grp,'-');
    
else
    
    % Just one group here
    if strcmp(str{val},'date')
        grp = statsDateConvert(meta.(str{val}));
    else
        grp = meta.(str{val});
    end
    
end



end

