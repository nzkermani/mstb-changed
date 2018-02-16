function [k, pval] = pvpair(pname,theVal,okargs,mfile,doPartialMatch)
%PVPAIR Helper function that looks for case insensitive partial/full
% matches of parameter names in a list of inputs and returns the
% parameter/value pair and matching number.
%
%   [K, PVAL] = PVPAIR(PNAME,THEVAL,OKARGS,CALLER) given input string
%   PNAME, and corresponding value, THEVAL, finds a partial matching name
%   in the OKARGS list. Returns K, the index of the match, and PVAL, the
%   parameter value. CALLER is the inserted string into the error
%   identifier, if any, i.e. bioinfo:CALLER:error_type. If pvpair is used
%   in a class method CALLER should be of the form class:method, otherwise
%   caller is the function name of the invoker.
%
%   [K, PVAL] = PVPAIR(PNAME,THEVAL,OKARGS,CALLER,false) finds a case
%   insensitive full match instead of the partial match. 

% Copyright 2008-2012 The MathWorks, Inc.

if nargin<5
    doPartialMatch = true;
end

if ~ischar(pname) || ~isrow(pname)
    msg = getString(message('bioinfo:pvpair:InvalidParameterName'));
    msgId = sprintf('bioinfo:%s:InvalidParameterName', mfile);
    x = MException(msgId,msg);
    x.throwAsCaller;
end

if doPartialMatch
    k = find(strncmpi(pname,okargs,numel(pname)));
else
    k = find(strcmpi(pname,okargs));
end

if numel(k) == 1
    pval = theVal;
    return
end

if isempty(k)
    msg = getString(message('bioinfo:pvpair:UnknownParameterName',pname));
    msgId = sprintf('bioinfo:%s:UnknownParameterName', mfile);
    x = MException(msgId,msg);
    x.throwAsCaller;

elseif length(k)>1
    msg = getString(message('bioinfo:pvpair:AmbiguousParameterName',pname));
    msgId = sprintf('bioinfo:%s:AmbiguousParameterName', mfile);
    x = MException(msgId,msg);
    x.throwAsCaller;
end
end %pvpair method
