function [dist,path,pred] = graphshortestpath(G,S,varargin)
%GRAPHSHORTESTPATH solves the shortest path problem in graph.
%
% [DIST,PATH,PRED] = GRAPHSHORTESTPATH(G,S) determines the single source
% shortest paths from node S to all other nodes in the graph G. Weights of
% the edges are all nonzero entries in the n-by-n adjacency matrix
% represented by the sparse matrix G. DIST are the n distances from source
% to every node (using Inf for non-reachable nodes and zero for the source
% node). The PATH contains the winning paths to every node, and PRED
% contains the predecessor nodes of the winning paths. 
% 
% [DIST,PATH,PRED] = GRAPHSHORTESTPATH(G,S,D) determines the single
% source-single destination shortest path from node S to node D.
% 
% GRAPHSHORTESTPATH(...,'METHOD',METHOD) selects the algorithm to use,
% options are:
%    'BFS'          - Breadth First Search, assumes all the weights are
%                     equal, edges are nonzero entries in the sparse matrix
%                     G. Time complexity is O(n+e).
%   ['Dijkstra']    - Assumes that weights of the edges are all positive
%                     values in the sparse matrix G. Time complexity is
%                     O(log(n)*e). 
%    'Bellman-Ford' - Assumes that weights of the edges are all nonzero
%                     entries in the sparse matrix G. Time complexity is
%                     O(n*e). 
%    'Acyclic'      - The input graph must be acyclic. Assumes that weights 
%                     of the edges are all nonzero entries in the sparse
%                     matrix G. Time complexity is O(n+e).
% 
% Note: n and e are number of nodes and edges respectively.
% 
% GRAPHSHORTESTPATH(...,'DIRECTED',false) indicates that the graph G is
% undirected, upper triangle of the sparse matrix is ignored. Default is
% true.
% 
% GRAPHSHORTESTPATH(...,'WEIGHTS',W) provides custom weights for the edges,
% useful to indicate zero valued weights. W is a column vector with one
% entry for every edge in G, traversed column-wise.
% 
% Examples:
%   % Create a directed graph with 6 nodes and 11 edges
%   W = [.41 .99 .51 .32 .15 .45 .38 .32 .36 .29 .21];
%   DG = sparse([6 1 2 2 3 4 4 5 5 6 1],[2 6 3 5 4 1 6 3 4 3 5],W)
%   h = view(biograph(DG,[],'ShowWeights','on'))
%   % Find shortest path from 1 to 6
%   [dist,path,pred] = graphshortestpath(DG,1,6)
%   % Mark the nodes and edges of the shortest path
%   set(h.Nodes(path),'Color',[1 0.4 0.4])
%   edges = getedgesbynodeid(h,get(h.Nodes(path),'ID'));
%   set(edges,'LineColor',[1 0 0])
%   set(edges,'LineWidth',1.5)
%
%   % Solving the previous problem for an undirected graph
%   UG = tril(DG + DG')
%   h = view(biograph(UG,[],'ShowArrows','off','ShowWeights','on'))
%   % Find the shortest path between node 1 and 6
%   [dist,path,pred] = graphshortestpath(UG,1,6,'directed',false)
%   % Mark the nodes and edges of the shortest path
%   set(h.Nodes(path),'Color',[1 0.4 0.4])
%   fowEdges = getedgesbynodeid(h,get(h.Nodes(path),'ID'));
%   revEdges = getedgesbynodeid(h,get(h.Nodes(fliplr(path)),'ID'));
%   edges = [fowEdges;revEdges];
%   set(edges,'LineColor',[1 0 0])
%   set(edges,'LineWidth',1.5)
%
% See also: GRAPHALLSHORTESTPATHS, GRAPHCONNCOMP, GRAPHISDAG,
% GRAPHISOMORPHISM, GRAPHISSPANTREE, GRAPHMAXFLOW, GRAPHMINSPANTREE,
% GRAPHPRED2PATH, GRAPHTOPOORDER, GRAPHTRAVERSE.
%
% References: 
%  [1]	E.W. Dijkstra "A note on two problems in connexion with graphs"
%       Numerische Mathematik, 1:269-271, 1959. 
%  [2]	R. Bellman "On a Routing Problem" Quarterly of Applied Mathematics,
%       16(1):87-90, 1958. 

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2007/10/05 18:32:12 $

algorithms = {'bfs','dijkstra','bellman-ford','acyclic'};
algorithmkeys = {'spb','dij','bel','spa'};
debug_level = 0;

% set defaults of optional input arguments
D = 1:length(G); % will return shortest path to all other nodes
W = []; % no custom weights
algorithm  = 2; % defaults to dijkstra
directed = true;

% find out signature of input arguments
if nargin>2 && isnumeric(varargin{1})
    D = varargin{1};
    varargin(1) = [];
end

% read in optional PV input arguments
nvarargin = numel(varargin);
if nvarargin
    if rem(nvarargin,2) == 1
        error('shortestpath:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'method','directed','weights'};
    for j=1:2:nvarargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error('shortestpath:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('shortestpath:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1 % 'method'
                    algorithm = strmatch(lower(pval),algorithms); %#ok
                    if isempty(algorithm) 
                        error('shortestpath:NotValidMethod',...
                              'String "%s" is not not a valid algorithm.',pval)
                    elseif numel(algorithm)>1
                         error('shortestpath:AmbiguousMethod',...
                              'String "%s" is ambiguous.',pval)
                    end
                case 2 % 'directed'
                    directed = opttf(pval,okargs{k},mfilename);
                case 3 % 'weights'
                    W = pval(:);
            end
        end
    end
end

% call the mex implementation of the graph algorithms
if nargout>1
    if isempty(W)
        [dist,pred] = graphalgs(algorithmkeys{algorithm},debug_level,directed,G,S);
    else
        [dist,pred] = graphalgs(algorithmkeys{algorithm},debug_level,directed,G,S,W);
    end    
else
    if isempty(W)
        dist = graphalgs(algorithmkeys{algorithm},debug_level,directed,G,S);
    else
        dist = graphalgs(algorithmkeys{algorithm},debug_level,directed,G,S,W);
    end
end

dist = dist(D);

% calculate paths
if nargout>1
    path = graphpred2path(pred,D);
end