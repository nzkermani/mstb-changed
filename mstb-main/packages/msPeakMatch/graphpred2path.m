function path = graphpred2path(pred,D)
%GRAPHPRED2PATH converts predecessor indices to paths.
%  
% PATH = GRAPHPRED2PATH(PRED,D) traces back a path by following the
% predecessor list in PRED starting at destination node D. PRED is a row
% vector of predecessor node indices and D is a scalar. The value of the 
% root (or source) node in PRED must be 0. PATH is a row vector listing the
% nodes from the root (or source) to D. If a NaN is found when following 
% the predecessor nodes, PRED2PATH returns an empty path. 
%
% When PRED and D are both row vectors, PATH is a row cell array with every
% column containing the path to the destination for every element in D. 
% When PRED is a matrix and D is a scalar, PATH is a column cell array with
% every row containing the path for every row in PRED.
% When PRED is a matrix and D is a row vector, PATH is a matrix cell array
% with every row containing the paths for the respective row in PRED, and
% every column containing the paths to the respective destination in D.
%
% If D is omitted, the paths to all the destinations are calculated for
% every predecessor listed in PRED.
%
% Example:
%   % Find the nodes from the root to one leaf in a phylogenetic tree, for
%   % instance, the GLR_HUMAN protein
%   
%   tr = phytreeread('pf00002.tree')
%   view(tr)
%   [CM,labels,dist] = getmatrix(tr);
%   root_loc = size(CM,1)
%   glr_loc = strmatch('GLR',labels)
%   [T,PRED] = graphminspantree(CM,root_loc);
%   PATH = graphpred2path(PRED,glr_loc)
%   labels(PATH)
%
% See also: GRAPHALLSHORTESTPATHS, GRAPHCONNCOMP, GRAPHISDAG,
% GRAPHISOMORPHISM, GRAPHISSPANTREE, GRAPHMAXFLOW, GRAPHMINSPANTREE,
% GRAPHSHORTESTPATH, GRAPHTOPOORDER, GRAPHTRAVERSE.  


%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2006/09/27 00:19:21 $

[numPreds,numNodes] = size(pred);
if nargin<2; D = 1:numNodes; end

% check input arguments
prednonan = ~isnan(pred);
if (any(pred(prednonan) < 0)) || (any(rem(pred(prednonan),1)))
    error('pred2path:InvalidIndices','Elements in PRED must be integers between 0 and the number of nodes (%d), or NaNs.',numNodes)
end
if (any(pred(prednonan) > numNodes)) 
    error('pred2path:InvalidIndicesTransposed','size(PRED) suggests there is (are) %d node(s), but there are elements in PRED larger, is PRED transposed?',numNodes)
end
if (any(D < 1)) || (any(D > numNodes)) || (any(rem(D,1)))
    error('pred2path:InvalidDestination','Elements in D must be integers between 1 and %d.',numNodes)
end

% intialize output
path = cell(numPreds,numel(D));
p = zeros(1,numNodes);

% trace back predecessors
for s = 1:numPreds
   for i = 1:numel(D)
       k = 1;
       n = D(i);
       p(k) = n;
       while (n~=0) && (k<=numNodes) && (~isnan(n))
           k = k+1;
           n = pred(s,n);
           p(k) = n;
       end
       if isnan(n)
           n = 0;
           path{s,i} = [];
       else
           path{s,i} = p(k-1:-1:1);
       end
       if n; break; end
   end
   if n; break; end
end
if n 
    error('pred2path:cycleDetected','PRED is invalid, potential cycle detected at node %d in row %d of PRED.',n,s)
end

if numel(path)==1; path = path{1}; end
