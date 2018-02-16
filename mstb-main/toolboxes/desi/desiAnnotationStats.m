function desiAnnotationStats(dpn)
%desiAnnotationStats - extract various parts from the guidata structure in
%order to provide a bit of info about intensities and stuff in the
%annotation files. This is designed, initially for QQQ data

% Annotations, and how many annotations are there?
[mask,histID,annoID] = desiAnnotationExtract(dpn,true);
numA = size(dpn.anno,1);

% Need to make an annotation figure with numbers on it to show which
% annotation is which
imageSubplot(dpn);

% Reshape the data
sz = size(dpn.d1.sp);
sp = reshape(dpn.d1.sp,[sz(1)*sz(2) sz(3)]);

% Somewhere to store the data
stMean = zeros(numA,size(sp,2));
stMedian = zeros(size(stMean));
annoNames = cell(1,numA);

% Loop through each annotation
for n = 1:numA
    
    fx = annoID == n;
    
    stMean(n,:) = nanmean(sp(fx,:),1);
    stMedian(n,:) = nanmedian(sp(fx,:),1);
    
    tmpName = ['Anno_' int2str(n) '_' dpn.anno{n,5}];
    spc = isstrprop(tmpName,'alphanum');
    tmpName(~spc) = '_';
    annoNames{1,n} = [tmpName '_Idx' int2str(n)];
    
end

% Make cell array of m/z values
mz = cell(size(sp,2),1);
for n = 1:size(mz,1)
    mz{n} = ['mz_' sprintf('%0.3f',dpn.d1.mz(n))];
    
    if n > 1 && strcmp(mz{n},mz{n-1})
        mz{n} = [mz{n} '*'];
    end            
end
    
% Make two tables...
tabMean = array2table(stMean',...
    'VariableNames',annoNames,...
    'RowNames',mz);    
tabMedian = array2table(stMedian',...
    'VariableNames',annoNames,...
    'RowNames',mz);
    

disp(['******** MEAN VALUES ********']);
disp(tabMean);
disp(['*****************************']);

disp(['******* MEDIAN VALUES *******']);
disp(tabMedian);
disp(['*****************************']);

% Return the tables to the workspace
assignin('base','tabMean',tabMean);
assignin('base','tabMedian',tabMedian);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imageSubplot(dpn)
% Make a good image...

tmp.fig = figure('Units','normalized',...
    'Position',[0.1 0.1 0.5 0.8]);
numA = size(dpn.anno,1);

tmp.ax1 = subplot(2,1,1);
imagesc(tmp.ax1,dpn.opt.coreg);
hold on;
for n = 1:numA
    [h] = redrawPatch(dpn.anno{n,6},...
        dpn.anno{n,7},dpn.anno{n,3},...
        [int2str(n) ': ' dpn.anno{n,5}]);
end
axis image
axis off

tmp.ax2 = subplot(2,1,2); imagesc(tmp.ax2,dpn.d1.img);
colormap(redbluecmap);
hold on;
for n = 1:numA
    [h] = redrawPatch(dpn.anno{n,8},...
        dpn.anno{n,9},dpn.anno{n,3},...
        [int2str(n) ': ' dpn.anno{n,5}]);
end
axis image
axis off

% Only link axes if the sizes are the same...
if size(dpn.d1.sp,1) == size(dpn.opt.coreg,1)
    linkaxes([tmp.ax1 tmp.ax2],'xy');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = redrawPatch(x,y,col,textLabel)
% Draw a patch in the CURRENT axes

if numel(x) == 1 && numel(y) == 1
    % Then this is a single scatter point of annotation
    h = scatter(x(1),y(1),200,col,'o','filled',...
        'MarkerEdgeColor','w');
    
else
    
    % Check that we have the annotations in the correct order for optical
    % and MS annotations
    x = [min(x) max(x) max(x) min(x)];
    y = [min(y) min(y) max(y) max(y)];

    % Then multiple points, so we patch them together
    x = x + [-0.5 0.5 0.5 -0.5];
    y = y + [-0.5 -0.5 0.5 0.5];
    
    h = patch(x,y,...
        col,...
        'EdgeColor','k',...
        'FaceColor',col,...
        'FaceAlpha',0.8,...
        'LineWidth',3);
end

% Add a text label
text(x(2)+5,mean(y(2:3)),textLabel,...
    'BackgroundColor',[0.7 0.7 0.7],...
    'FontSize',16,...
    'FontWeight','bold',...
    'Clipping','on');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%