function [ output_args ] = mzDeviation(all)
%mzDeviation - given the cell output from annotation2text, produce a
%diagnostics style plot to show how the mz values deviate when provided
%with a annotations, i.e. plot the true and experimentally measured values
%for inspection of e.g. systematic drift.

% First we should import the Ottmar-approved list of calibration ions from
% the Excel spreadsheet
omgPath = ['packages' filesep 'msDatabases' filesep 'mz4calib.xlsx'];

[calib.ions,calib.name,~] = xlsread(omgPath,'Sheet1');

% Best to convert to a column vector with numbers rather than the current
% text based cell format.
[op] = convert2vector(all,'notfull');

% Find the calibration ions in the assignment list..
[inds,cname] = findCalibrationIons(op,calib);


% Calculate the ppm deviation
ppm = 1e6 * (op(:,1) - op(:,2)) ./ op(:,2);

% Now the scatter using just the omg ions
figure('Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5],...
    'Toolbar','figure',...
    'Menubar','none');
hold on;

scatter(op(inds(:,1),2),ppm(inds(:,1),1),80,'blue','o','filled');
for n = 1:size(inds,1)
    text(op(inds(n,1),2),ppm(inds(n,1),1),['  --- ' cname{1,n}],'FontSize',12);
end
xlabel('True m/z', 'FontSize',16,'FontWeight','bold');
ylabel('m/z deviation / ppm', 'FontSize',16,'FontWeight','bold');
title('Key Diagnostic Ions','FontSize',18);
box on;




return

figure; hold on;
scatter(op(:,2),op(:,1),50,'blue','o','filled');
line([600 1000],[600 1000],'Color','black');

figure; hold on;
scatter(op(:,2),ppm,50,'blue','o','filled');


% Do LOESS smoothing...
sm1 = smooth(op(:,2),ppm,100,'rlowess');
sm2 = smooth(op(:,2),ppm,100,'rloess');

plot(op(:,2),sm1,'red','LineWidth',4);
plot(op(:,2),sm2,'cyan','LineWidth',4);



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x] = convert2vector(all,mode)
% Convert the cell into a set of numbers for scatter plotting. You can
% chose to include instances of multiple assignments or not depending on
% the 'mode' (explain later)

% How many variables?
numV = size(all,1);

% Create empty matrix that is large enough to hold it all...
x = zeros(numV,2);

% Now run through the list adding into x
i = 0;
for n = 1:numV
    
    % Check that there is actually an assignment for this exp mz value
    if ~isempty(all{n,2})
        
        % Now see if there is more than one true mz value (i.e. multiple
        % assignments)
        cm = strfind([all{n,2} ','],',');
        
        if numel(cm) == 1
            
            i = i + 1;            
            x(i,:) = [all{n,1} str2double(all{n,2})];
            
        elseif numel(cm) > 1 && strcmp(mode,'full')
            % So here there are multiple assignments, which you may wish to
            % include or not...            
            lst = 1;
            for r = 1:numel(cm)
                i = i + 1;
                x(i,1) = all{n,1};
                x(i,2) = str2double(all{n,2}(lst:cm(r)-1));                
                lst = cm(r)+1;                
            end
            
        else
            % Do nothing - this will happen when mode ~= 'full'        
        end
    end    
end
       
% Trim the full list
x = x(1:i,:);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inds,name] = findCalibrationIons(op,calib)
% Find the calibration ions in the annotation list thing

numC = numel(calib.ions);
tol = 0.001;
inds = zeros(numC,1);

for n = 1:numC
    
    % Find true values that are v. v. close...
    [fx,~] = find(op(:,2)-tol < calib.ions(n) & op(:,2)+tol > calib.ions(n));
    
    % Include only if an unambiguous match
    if numel(fx) == 1
        inds(n,1) = fx;
    end    
end

[fx,~] = find(inds ~= 0);

inds = [inds(fx,1) calib.ions(fx)'];
name = calib.name(1,fx);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%