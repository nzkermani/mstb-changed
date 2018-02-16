function [ output_args ] = mtspAnnotationCount(res)
%mtspAnnotationCount - determine the quantity and quality of annotations
%from the online engine.

% It is clear that some files have so few (good) annotations that it would
% be an interesting side analysis, especially if we compare to the
% annotations of the lipid database that we maintain at Imperial.


numF = size(res,1);
sz = zeros(numF,2);
for n = 1:numF
   
    % File's annotations
    tmp = res{n};
    
    % Size
    sz(n,1) = size(tmp,1);
    
    % Quantity of FDR passes (0.1 threshold default)
    fx = strcmpi(tmp(:,10),'true');
    sz(n,2) = sum(fx);
    
    
end

lim = max(sz(:))

figure; 
subplot(1,2,1); hist(sz(:,1),0:10:lim);
xlim([-10 lim+10]);
ylim([-0.1 15]);
set(gca,'FontSize',14);
title('Annotations: all','FontSize',20,'FontWeight','bold');
xlabel('Number of Annotations','FontWeight','bold','FontSize',18);
ylabel('Quantity of Files','FontWeight','bold','FontSize',18);

subplot(1,2,2); hist(sz(:,2),0:10:lim);
xlim([-10 lim+10]);
ylim([-0.1 15]);
set(gca,'FontSize',14);
title('Annotations: FDR < 0.1','FontSize',20,'FontWeight','bold');
xlabel('Number of Annotations','FontWeight','bold','FontSize',18);
ylabel('Quantity of Files','FontWeight','bold','FontSize',18);


end

