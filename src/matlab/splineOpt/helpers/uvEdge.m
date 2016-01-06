function [ splines, eindices ] = uvEdge( uInd, vInd, V, E, adjList )
%UVEDGE finds the spline edge with vertices u and v, if it exists.
%
% INPUTS
% uInd: (integer) index of the vertex endpoint u of interest.
% vInd: (integer) index of the vertex endpoint v of interest.
% V: (cell array of (dx1) vertex coordinates) vertex set.
% E: (cell array of spline structs) of splines.
% adjList: (Nx1 cell array): for index i, contains a (2xK) matrix of K
% vertex-edge index pairs of vertices and edges adjacent to vertex i.
%
% OUTPUTS
% splines: (cell of struct) splines connected by u and v, or {} if such a
% spline doesn't exist.
% eindices: (integer array) indices of splines with endpoints u and v, or
% [] if such a spline doesn't exist.
%
% @author Roger Zou
% @date 6/13/15

% get neighboring edges of both u and v
[edgeIa] = vNeighborE(uInd, V, E, adjList);
[edgeIb] = vNeighborE(vInd, V, E, adjList);

% find the intersection and output
eindices = intersect(edgeIa, edgeIb);
if isempty(eindices)
    splines = {};
else
    splines = E(eindices);
end

end
