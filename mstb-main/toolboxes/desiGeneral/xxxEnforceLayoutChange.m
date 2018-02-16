function xxxEnforceLayoutChange(fig,mode,option)
%xxxEnforceLayoutChange

% Need to provide a good input
if ~strcmp(option,'on') && ~strcmp(option,'off')
    return
end

% Determine current state
cus = get(fig.tb.layout,'State');

% Then the function will do nothing as it is as requested already
if strcmp(cus,option)    
    return
end

% Otherwise, we switch...
if strcmp(cus,'on')
    set(fig.tb.layout,'State','off');
        

else
    set(fig.tb.layout,'State','on');
    
end

% Then run either the desi/desiPosNeg function to change the layout
if strcmp(mode,'single')
    desiChangeLayout(fig.tb.layout,[],fig);
else
    dpnChangeLayout(fig.tb.layout,[],fig);
end

end

