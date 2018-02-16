function rawH5plotMeanSpectra(spec,resamp,picked)
%rawH5plotMeanSpectra - blah blah blah

numS = size(spec,1);
cols = parula(numS);

figure; hold on;

for n = 1:numS
    
%     plot(spec{n}(:,1),spec{n}(:,2),...
%         'Color',cols(n,:),...
%         'LineWidth',2);
%     
%     plot(resamp.reMZ,resamp.reSp(n,:),...
%         'Color',cols(n,:),...
%         'LineWidth',2,...
%         'LineStyle','-');
    
    plot(picked.refSpec(:,1),picked.al(:,n),...
        'Color',cols(n,:),...
        'LineWidth',2,...
        'LineStyle','-');
    
    

end



end

