function [y,scaleFac] = jsmNormalise(x,method,thresh,rescale,grouping,normRef)
% x - sparse matrix of data
% method - one of the named methods
% thresh - a threshold below which data don't count
% rescale - currently not used
% grouping - necessary for the slide/section methods, differential
% scaling factors according to the 
%tt = tic;

% Determine the original mean intensity of the raw data, which will be used
% for re-scaling purposes
if rescale
    rescaleValue = nanmax(nanmean(x,1))
end

% Check that is of the correct size (i.e. number of variables)
if nargin > 5
    if numel(normRef) ~= size(x,2)
        error('normRef is incorrectly sized');
    %else
        %normRef = [];
        %error('Not a good normRef');
    end
else
    normRef = [];
end

% Some methods run by groupings, e.g. s-tic, s-med. Check that we have them
% if needed, and then calculate the unique groups within.
switch lower(method)
    case {'s-tic','s-med'}
        if isempty(grouping)
            grouping = ones(size(x,1),1);
            warning('stic method won''t work!');
        end
        
        % Determine unique groupings and indices
        [unq,~,unqI] = unique(grouping);        
        
    otherwise
        % Grouping not required
end

% A threshold?
if isempty(thresh)
    thresh = 0;
end

% First we need to determine the scaling factor
switch lower(method)
    case 'tic'
        % Find the sum of each pixel / spectrum
        scaleFac = nansum(x,2);

    case 's-tic'
        % This creates a normalisation factor based on the sum of the
        % entire tissue section rather than individual pixels. So we need
        % to tell it which pixels belong to which files...              
        scaleFac = NaN(size(x,1),1);
        
        % Run through each level to get the total sum over the grouping,
        % which is typically the tissue section...
        for n = 1:numel(unq)            
            fx = unqI == n;
            tmp = full(nansum(nansum(x(fx,:))));
            scaleFac(fx,1) = tmp;
        end
        
    case {'pqn-median','pqn-med','pqn-mean','pqn-mean-all'}
        % A slightly different implementation of the PQN policy
        tmpRef = lower(method(5:end));
        
        % Find values above the threshold
        if issparse(x)
            tmp = sparse(NaN(size(x)));
        else
            tmp = NaN(size(x));
        end
        tmp(x > thresh) = x(x > thresh);

        % Determine the reference spectrum, unless it is provided
        if isempty(normRef)
            switch tmpRef
                case 'mean'
                    normRef = nanmean(tmp,1);
                otherwise
                    normRef = nanmedian(tmp,1);
            end
        end
        
        % Now determine the fold changes
        fldchn = bsxfun(@rdivide,tmp,normRef);
        
        % Now calculate the median of these fold changes for each sample
        scaleFac = nanmedian(fldchn,2);

    case 'uve'
        % Find the variables with the largest mean intensity
        mnn         = nanmedian(x,1);
        [~,bpk]     = max(mnn);
        scaleFac    = x(:,bpk);

    case 'med'
        % Find the values above the threshold
        if issparse(x)
            tmp = sparse(NaN(size(x)));
        else
            tmp = NaN(size(x));
        end        
        tmp(x > thresh) = x(x > thresh);
        
        scaleFac = nanmedian(tmp,2);
        
    case 's-med'
        % This calculates a median value as per the grouping variable
        % provided, akin to finding the sum over the whole section in stic
        scaleFac = NaN(size(x,1),1);
        
        % Run through each level to get the total sum over the grouping,
        % which is typically the tissue section...
        for n = 1:numel(unq)            
            fx = unqI == n;
            tmp = x(fx,:);        
            tmp = median(tmp(tmp > thresh));
            scaleFac(fx,1) = full(tmp);
        end
        
    case 'med-fc'
        % Calculate the median over the entire set of spectra
        if issparse(x)
            tmp = sparse(NaN(size(x)));
        else
            tmp = NaN(size(x));
        end        
        tmp(x > thresh) = x(x > thresh);
        
        scaleFac = nanmedian(tmp,1);

    case 'mean-fc'
        % Calculate the mean over the entire set of spectra
        if issparse(x)
            tmp = sparse(NaN(size(x)));
        else
            tmp = NaN(size(x));
        end        
        tmp(x > thresh) = x(x > thresh);
        
        scaleFac = nanmean(tmp,1);

    case 'raw'
        % Do nothing...
        y = x;
        scaleFac = [];
        disp('no norm');
        return

    otherwise
        y = x;
        scaleFac = [];
        disp('no norm');
        return
end
    
% Set all Inf, Nan, 0 to values of 1
scaleFac(isinf(scaleFac))   = 1;
scaleFac(isnan(scaleFac))   = 1;
scaleFac(scaleFac == 0)     = 1;

% Are there any variables that sum to 0, i.e. fully empty?
zzz = sum(x,1) == 0;

% Now that we have an appropriate scaling factor, we use it!
y = bsxfun(@rdivide,x,scaleFac);

% Set empty variables to zero
y(:,zzz) = 0;

% Consider performing rescaling in order to put the spectra back to an
% intensity-appropriate level, rather than being stuck around 0-1
if rescale
    y = y * rescaleValue;
end

%toc(tt);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
