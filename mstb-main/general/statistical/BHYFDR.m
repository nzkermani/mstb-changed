function [pThrICPC,pThrNC] = BHYFDR(pvalues,aplha)

%% The function is used to find the p-value threshold to control the given alpha false 
%% discovery rate under dependence assumptions of metabolite peaks via the Benjamini–Hochberg–Yekutieli
%% procedure
%%  Input: 
%           pvalues - pvalues for testing of no singificance of metabolite
%                   changes between treatment groups
%           alpha   - false discovery rate
%%  Output: 
%           pThrICPC    - p threshold calculated by assuming that metabolite
%                         variables are independent or positively correlated
%           pThrNC      - p threshold calculated by assuming an arbitrary
%                         correlation structure between metabolite
%                         variables
%% Author: Kirill Veselkov, Imperial College London, 2010

dims     = size(pvalues);
dimIndcs = (dims>1);
if sum(dimIndcs)==2
    if dims(1)>dims(2)
        pvalues = pvalues';
    end
    nMthds   = size(pvalues,1);
    pThrICPC = zeros(1,nMthds);
    pThrNC   = zeros(1,nMthds);
    for iMthd = 1:nMthds
        [pThrICPC(iMthd),pThrNC(iMthd)] = BHYFDR(pvalues(iMthd,:),aplha);
    end
    return;
end

nVrbls             = length(pvalues);
orderedPvls        = sort(pvalues);
k                  = 1:nVrbls;

% find the largest k such that
cumsumVrbls        = sum(1./(1:nVrbls));
imax = find(orderedPvls<=aplha.*k./(nVrbls.*cumsumVrbls),1,'last');
if ~isempty(imax)
    pThrNC = orderedPvls(imax);
else
    pThrNC = 0;
end

imax = find(orderedPvls<=aplha.*k./nVrbls,1,'last');
if ~isempty(imax)
    pThrICPC = orderedPvls(imax);
else
    pThrICPC = 0;
end

return;