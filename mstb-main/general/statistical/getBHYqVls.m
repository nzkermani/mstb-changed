function [qVlsICPC,qVlsNC] = getBHYqVls(pvalues,alpha)
%% getBHYqVls: Calculation of the p-value threshold to control the given alpha false 
%% discovery rate under dependence assumptions via the Benjamini–Hochberg–Yekutieli
%% procedure
%%  Input: 
%           pvalues [methods x variables] - pvalues for testing of no singificance of metabolite
%                                           changes between treatment groups
%           alpha   - false discovery rate
%%  Output: 
%           qVlsICPC    - FDR values assuming that metabolite variables are independent or positively correlated
%           qVlsNC      - FDR values assuming that metabolite variables are arbitrary correlated
%% Author: Kirill Veselkov, Imperial College London, 2010

dims     = size(pvalues);
dimIndcs = (dims>1);
if sum(dimIndcs)==2
    if dims(1)>dims(2)
        pvalues = pvalues';
    end
    [nMthds,nVrbls]   = size(pvalues);
    qVlsICPC = zeros(nMthds,nVrbls);
    qVlsNC   = zeros(nMthds,nVrbls);
    for iMthd = 1:nMthds
        [qVlsICPC(iMthd,:),qVlsNC(iMthd,:)] = getBHYqVls(pvalues(iMthd,:),alpha);
    end
    return;
end


nVrbls                   = length(pvalues);
[sortPvls,sortIndcs]     = sort(pvalues);
k                        = 1:nVrbls;
cumsumVrbls              = sum(1./(1:nVrbls));

qVlsNC                   = NaN(1,nVrbls);
qVlsNCsorted             = sortPvls.*nVrbls.*cumsumVrbls./k;
qVlsNC(sortIndcs)        = qVlsNCsorted;
qVlsNC(qVlsNC>alpha)     = 1;

qVlsICPC                 = NaN(1,nVrbls);
qVlsICPCsorted           = sortPvls.*nVrbls./k;
qVlsICPC(sortIndcs)      = qVlsICPCsorted;
qVlsICPC(qVlsICPC>alpha) = 1;
return;