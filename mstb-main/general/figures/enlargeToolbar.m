function enlargeToolbar(hTB, iconsizeTB, fontsize)
%% Enlarge toolbar

if nargin<2
    iconsizeTB = 60;
end

if nargin<3
    fontsize = 12;
end
drawnow;

jToolbar   = get(get(hTB ,'JavaContainer'),'ComponentPeer');
jToolbar.setPreferredSize(java.awt.Dimension(iconsizeTB*1.3,iconsizeTB*1.3));
    
newSize    = java.awt.Dimension(iconsizeTB*1.3,iconsizeTB*1.3);
numButtons = jToolbar.getComponentCount;
for i = 1:numButtons
    toolTipString = jToolbar.getComponent(i-1).getToolTipText;
    %toolTipString = string(toolTipString);
    toolTipString = char(toolTipString);
    if ~isempty(toolTipString)
        indcs = strfind(toolTipString,':');
        if ~isempty(indcs)
            jToolbar.getComponent(i-1).setToolTipText(toolTipString(1:indcs(end)-1));
            jToolbar.getComponent(i-1).setText(toolTipString(indcs(end)+1:end));
            jToolbar.getComponent(i-1).setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
            jToolbar.getComponent(i-1).setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
            jToolbar.getComponent(i-1).setFont(java.awt.Font('Garamond',java.awt.Font.PLAIN,fontsize));
        end
    end
    jToolbar.getComponent(i-1).setMaximumSize(newSize);
end
jToolbar.revalidate; % refresh/update the displayed toolbar clear jToolbar; drawnow; set(msimage.mainfig.hTB,'Visible','on');
clear jToolbar; 
set(hTB,'Visible','on');
end