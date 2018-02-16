function [ tab ] = mtspCSVcompare(res,dpn)
%mtspCSVcompare - compare the annotations coming from two annotations of th
%same file.  Expect there to be missing annotations but mostly similar.
%Want to compare actual annotations quantites, scores and other things like
%that. Use the chemical formulae + adduct to match as these will be unique.

% Originally designed to compare only 2 files!
old = 'ICL//LNTO41_17_5';
new = 'LNTO41_17_5-RECAL';

% Generate the matched table
[nh,tab] = makeTable(res,old,new);

% Find the ions only in the original file
%ionExtractImageMissing(tab,dpn.mtspOLD)
ionExtractImageGained(tab,dpn.mtsp)

return

% Draw the MSM graph
msmGraph(tab);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ionExtractImageMissing(tab,data)

% Find the ions unique to the original file
fdrs = cell2mat(tab(:,7:8))

fdrs(isnan(fdrs)) = Inf;

fx = fdrs(:,1) <= 0.1;
fy = fdrs(:,2) > 0.1;

% Temporary table and noew image storage
newtab = tab(fx&fy,:);
allI = cell(size(newtab,1),1);

% Can we match the data...
for n = 1:size(newtab,1)
    
    thisForm = newtab{n,1};
    thisAdct = newtab{n,2};
    
    % Match to existing formulae
    fx = strcmp(data.annos(:,1),thisForm);
    fy = strcmp(data.annos(:,2),thisAdct);
    fz = fx & fy;
    if sum(fz) == 1
        allI{n,1} = data.sp(:,:,fz);        
    else
        allI{n,1} = zeros(size(data.sp,1),size(data.sp,2));
        disp('not found');
    end
end

imzmlIonTile(allI,1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ionExtractImageGained(tab,data)

% Find the ions unique to the original file
fdrs = cell2mat(tab(:,7:8))

fdrs(isnan(fdrs)) = Inf;

fx = fdrs(:,2) <= 0.1;
fy = fdrs(:,1) > 0.1;

% Temporary table and noew image storage
newtab = tab(fx&fy,:);
allI = cell(size(newtab,1),1);

% Can we match the data...
for n = 1:size(newtab,1)
    
    thisForm = newtab{n,1};
    thisAdct = newtab{n,2};
    
    % Match to existing formulae
    fx = strcmp(data.annos(:,1),thisForm);
    fy = strcmp(data.annos(:,2),thisAdct);
    fz = fx & fy;
    if sum(fz) == 1
        allI{n,1} = data.sp(:,:,fz);        
    else
        allI{n,1} = zeros(size(data.sp,1),size(data.sp,2));
        disp('not found');
    end
end

%imzmlIonTile(allI,1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msmGraph(tab)

% Find common annotations with fdr <= 0.1
fdrs = cell2mat(tab(:,7:8));
fdrs(isnan(fdrs)) = Inf;

% FInd annotations with both FDRs <= 0.1
fx = max(fdrs,[],2) <= 0.1;

msms = cell2mat(tab(fx,5:6));
ff = msms(:,1) >= msms(:,2);


figure; hold on;
scatter(msms(ff,1),msms(ff,2),80,'red','o','filled','MarkerEdgeColor','k');
scatter(msms(~ff,1),msms(~ff,2),80,'green','o','filled','MarkerEdgeColor','k');
line([0 1],[0 1],'LineStyle','--');
axis equal
xlim([0 1]);
ylim([0 1]);
box on;
set(gca,'FontSize',14);
xlabel('MSM (original)','FontWeight','bold');
ylabel('MSM (recalibrated)','FontWeight','bold');



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compareTo(tab)

% Compare to 0.1 fdr annotations from the old list
fx = cell2mat(tab(:,7)) <= 0.1;
%compareTo(tab(fx,:));

fy = isnan(cell2mat(tab(:,3)));
fz = cell2mat(tab(:,8)) <= 0.1;
compareTo(tab(fy&fz,:));

msms = cell2mat(ip(:,5:6));
fdrs = cell2mat(ip(:,7:8));
fdrs(isnan(fdrs)) = Inf;

fx = fdrs(:,2) <= 0.1;

ip(fx,:)
sum(fx)

ff = cell2mat(ip(fx,5));

[min(ff) median(ff) mean(ff) max(ff)]

figure; boxplot(cell2mat(ip(fx,5)))

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nh,tab] = makeTable(res,old,new)

% Find unique annotations...
allA = classMany2One([res(:,4) res(:,5)]);
[unq,~,ind] = unique(allA);
numA = numel(unq);

% New table to contain all unique annotations - larger than required but
% will be trimmed down
nh = {'Formula','Adduct',...
    'mzOld','mzNew',...
    'msmOld','msmNew',...
    'fdrOld','fdrNew',...
    'rhoOld','rhoNew'};
tab = cell(numA,numel(nh));

% Now let's put the ones in a table that match...
for n = 1:numA
    
    % This annotation
    fx = ind == n;    
    tmp = res(fx,:);
    
    % Is there an 'old' entry?
    oe = strcmp(tmp(:,2),old);
    if any(oe)        
        tab(n,1) = tmp(oe,4);
        tab(n,2) = tmp(oe,5);
        tab(n,[3 5 7 9]) = [tmp(oe,6) tmp(oe,7) tmp(oe,8) {tmp(oe,[9:11])}];
    else
        tab(n,[3 5 7 9]) = {NaN NaN NaN NaN};
    end
    
    % Is there an old entry?
    ne = strcmp(tmp(:,2),new);
    if any(ne)
        tab(n,1) = tmp(ne,4);
        tab(n,2) = tmp(ne,5);
        tab(n,[4 6 8 10]) = [tmp(ne,6) tmp(ne,7) tmp(ne,8) {tmp(ne,[9:11])}];
    else
        tab(n,[4 6 8 10]) = {NaN NaN NaN NaN};
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


