function [ edgeInds, edges, vertices ] = vNeighborE( vInd, V, E, adjList )
%VNEIGHBORE For a vertex, find adjacent edges (and vertices).
%
% INPUTS
% vInd: (integer) index of vertex in V.
% V: (cell array) of vertex coordinates
% E: (cell of spline structs) of splines, some of which may be
% adjacent to vertex v
% adjList: (Nx1 cell array): for index i, contains a (2xK) matrix of K
% vertex-edge index pairs of vertices and edges adjacent to vertex i.
%
% OUTPUTS
% edgeInds: (Kx1 integer array) edge indices of edges adjacent to v.
% edges: (cell of spline structs) splines adjacent to v.
% vertices: (cell of (d x 1) vectors) of adjacent vertices to v.
%
% @author Roger Zou
% @date 6/13/15

adjlisti = adjList{vInd};
edgeInds = adjlisti(2,:)';
if nargout > 1
    edges = E(edgeInds);
end
if nargout > 2
    vertices = V( adjlisti(1,:) );
end
