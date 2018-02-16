function dpnIonImage(~,~,ax,sp)
% Change the ion image displayed in the MS1/MS2 axes

% Determine what kind of image is provided.  If multidimensional then we
% just sum over the dimeinsions.
if numel(size(sp)) == 3
    if size(sp,3) == 3
        % Use as is
        sp = imScale(sp);
    else
        sp = nansum(sp,3);
    end
end

% Quick check to decide if we have been provided with axes/image handles,
% or just the axes.
if numel(ax) == 1
    
    % Get the children from the axes parent
    allCh = get(ax(1),'Children');
    
    % If there are absolutley none, then just make one
    if numel(allCh) == 0
        
        ax(2) = imagesc(ones(1,1,3),'Parent',ax(1));
    
    else
        
        % Find the type to see if there are any images there
        typCh = get(allCh,'Type');
        fx = strcmp(typCh,'image');
        
        % These are the ones to be deleted once we've finished
        f0 = allCh(~fx);
        
        % Find the image(s)
        fx = find(fx == 1);
        
        % If one, then simply use that. If none then make one. If more than
        % one, then what?
        if numel(fx) == 1     
            ax(2) = allCh(fx);
        elseif numel(fx) == 0
            ax(2) = imagesc(ones(1,1,3),'Parent',ax(1));
        else
            error('don''t have a handle');            
        end
        
        % Delete the extra children
        delete(f0);
    end
    
end

% The second element refers to the image, not the axes
set(ax(2),'CData',sp,...
    'Visible','on',...
    'Clipping','on');

% Default colormap...
colormap(redbluecmap);
% try
%     caxis([min(sp(:)) max(sp(:))]);
% catch err
%     err
% end

% This changes the axes
set(ax(1),'XTick',[],...
    'YTick',[],...    
    'XLim',[1 size(sp,2)],...    
    'YLim',[1 size(sp,1)],...
    'Visible','on',...
    'YDir','reverse',...
    'LineWidth',5);

set(ax(1),'LineWidth',5);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%