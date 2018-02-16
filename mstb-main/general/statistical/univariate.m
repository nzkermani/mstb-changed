function [ pq ] = univariate(data,grp,varargin)
%univariate - perform univariate analysis and q-value correction on a data
%matrix of multiple variables.

% Read in the options
[opts] = readArgsData(varargin,grp);

% Create a matrix to store the p values
numV = size(data,2);
pval = ones(1,numV);

% Define the function that we will use...
switch opts.test
    case 'ttest'
        fh = @uvttest;
    case 'ttest2'
        fh = @uvttest2;
    case 'anova'
        fh = @uvanova1;
    case {'kw','kruskal','kruskalwallis','kruskal-wallis'};
        fh = @uvkw;
    otherwise
        error('Unknown request');
end

% Loop through and perform
wb = waitbar(0,'UNIVARIATE');
warning off all
for n = 1:numV    
    pval(n) = fh(data(:,n),grp,opts.alpha);
    waitbar(n/numV,wb);
end
delete(wb);
warning on all

% Perform Q value correction
qval = getBHYqVls(pval,opts.alpha);

% Export the p and q values
pq = [pval' qval'];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = readArgsData(argsin,grp)
% Read the arguments and then the data if it wasn't passed

% Determine the number of groups...
[unq,~,~] = unique(grp);
numG = numel(unq);

% Define the defaults here
opts.alpha = 0.05;
switch numG
    case 1
        opts.test = 'ttest';
    case 2
        opts.test = 'ttest2';
    otherwise
        opts.test = 'anova';
end
opts.qValMethod = 'bhy';


% Run through each pair
numA = length(argsin);
for n = 1:2:numA
    
    % This is the 'name'
    name = lower(argsin{n});
    
    % Switch through the various options
    switch name
        
        case {'test','method'}
            opts.test = lower(argsin{n+1});
            
        case 'alpha'
            opts.alpha = argsin{n+1};
            
        case {'pq','pqcorr','qval'}
            opts.qValMethod = argsin{n+1};
            
        otherwise
            % There is no otherwise
            disp(['Unknown option: ' name]);
               
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p] = uvttest(y,~,alpha)
% Write this eventually

p = ttest(y,'Alpha',alpha);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p] = uvttest2(y,grp,alpha)
% Run the ttest2

[~,~,unq] = unique(grp);

[~,p] = ttest2(y(unq),y(~unq),'Alpha',alpha);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p] = uvanova1(y,grp,~)

p = anova1(y,grp,'off');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p] = uvkw(y,grp,~)

p = kruskalwallis(y,grp,'off');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


