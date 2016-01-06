function [ pathCell ] = embryoGTGraph( G )
%GROUNDTRUTHGRAPH Construct graph embryo model from ground truth.
%
% INPUTS
% G: (logical matrix) ground truth image of embryo with M edges.
%
% OUTPUTS
% pathCell: (M x 1) cell array of edge samplings. Entry i, a sampling of
% the i-th edge, is a (2 x n) matrix of n samples in R^2.
%
% @author Roger Zou
% @date 8/13/15

% Find all vertices and edges of graph with Grid object from ground truth G
truthGrid = javaObjectEDT('TruthGrid', G);
allEdges = truthGrid.getAllEdges();
srcV = squeeze(allEdges(1,:,:));
dstV = squeeze(allEdges(2,:,:));
M = size(srcV,1);

% get source, dest vertices of unique edges (remove duplicates)
allptsL = [srcV, dstV];
allptsR = [dstV, srcV];
allptsnorm = [ sqrt(sum(srcV.^2,2)), sqrt(sum(dstV.^2,2)) ];
[~, I] = sort(allptsnorm, 2);
sortedpts = NaN(M, 4);
orderL = I(:,1)==ones(M,1);
orderR = I(:,2)==ones(M,1);
sortedpts(orderL,:) = allptsL(orderL,:);
sortedpts(orderR,:) = allptsR(orderR,:);
[~, ia, ~] = unique(sortedpts, 'rows');
srcV = srcV(ia,:);
dstV = dstV(ia,:);

% obtain path sampling from ground truth
xE = ones(size(G)) * sqrt(2);
tE = ones(size(G));
% use pixel distance for edge costs
edgeCosts = cat(3, tE, xE, tE, xE);
% very large vertex costs for vertices, 1 for edge, 0 otherwise.
vertCosts = double(imcomplement(G));
srcVind = sub2ind(size(G), srcV(:,1)+1, srcV(:,2)+1);
dstVind = sub2ind(size(G), dstV(:,1)+1, dstV(:,2)+1);
Vind = union(srcVind, dstVind);
vertCosts(Vind) = sum(vertCosts(:))*5;
% construct ImageGraph object
iGraph = javaObjectEDT('ImageGraph', vertCosts, edgeCosts, G);

% retrieve list of samples of all edges
pathList = iGraph.getEdgePaths(int32(srcV), int32(dstV));
pathCell = cell(pathList.size(), 1);
for i = 1:pathList.size()
    % get path sample of each ground truth edge curve
    path = pathList.get(i-1) + 1;
    pathCell{i} = path';
end

end
