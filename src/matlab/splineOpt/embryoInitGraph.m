function [ V, E, adjList, faceList, Gamma ] = embryoInitGraph( GT, spacing, optInitSpline, parallel )
%EMBRYOINITGRAPH Compute initial splines and vertices from ground truth.
%
% INPUTS
% GT: ({0,1}^2 matrix) Ground truth of embryo image.
% spacing: (scalar) spacing between each control point in pixels.
% optInitSpline: function handle corresponding to edge type.
% parallel: parallelization toggle.
%
% OUTPUTS
% V: (cell array of (x,y) coord) of all junctions in GT, aka the vertex set.
% E: (cell array of struct) of all splines in GT, aka the edge set.
% adjList: (Nx1 cell array): for index i, contains a (2xK) matrix of K
% vertex-edge index pairs of vertices and edges adjacent to vertex i.
% faceList:
% Gamma: (cell array of R^2 matrix) of the original samplings for each
% corresponding spline in Splines.
%
% @author Roger Zou
% @date 8/13/15

% read ground truth and retrieve sampling of all edges
[pathCell] = embryoGTGraph(GT);

% get dimensions
M = numel(pathCell);

%% iterate over each edge, compute optimal spline.
E = cell(M, 1);
Gamma = cell(M, 1);
if parallel
    parfor uInd=1:M
        % get edge sampling gamma
        gamma = pathCell{uInd};
        gamma = [gamma(2,:); gamma(1,:)];
        Gamma{uInd} = gamma;
        % dynamically compute number of control points
        k = nCtrlPts(gamma, spacing);
        % create optimal spline from sampling
        s = optInitSpline(gamma, 3, k);
        E{uInd} = s;
    end
else
    for uInd=1:M
        % get edge sampling gamma
        gamma = pathCell{uInd};
        gamma = [gamma(2,:); gamma(1,:)];
        Gamma{uInd} = gamma;
        % dynamically compute number of control points
        k = nCtrlPts(gamma, spacing);
        % create optimal spline from sampling
        s = optInitSpline(gamma, 3, k);
        E{uInd} = s;
    end
end

srcV = zeros(M,2);
dstV = zeros(M,2);
for ii=1:M
    srcV(ii,:) = fliplr(E{ii}.control(:,1)');
    dstV(ii,:) = fliplr(E{ii}.control(:,end)');
end

% Construct Graph object from vertex/edge matrices and obtain vertex info
curGraph = javaObjectEDT('JunctionGraph', srcV, dstV);
Varray = double(curGraph.V());
Varray = [Varray(:,2), Varray(:,1)];
V = num2cell(Varray', 1)';

% construct adjacency list if required
if nargout > 2
    javaAdjList = curGraph.getAdjList(srcV, dstV);
    adjList = cell(javaAdjList.size(),1);
    for ii=1:javaAdjList.size()
        adjList{ii} = javaAdjList.get(ii-1)';
    end
end

% construct face list if required
if nargout > 3
    javaFaceList = curGraph.getFaceList(srcV, dstV);
    faceList = cell(javaFaceList.size(),1);
    for ii=1:javaFaceList.size()
        faceList{ii} = javaFaceList.get(ii-1)';
    end
end

end
