function [ g ] = DgcostDv( uInd, vInd, structA, structB, structC, structD )
%DGCOSTDV Gradient of the graph cost wrt to vertex v.
%
% INPUTS
% uInd: (integer) index of the endpoint of interest in image I1.
% vInd: (integer) index of the endpoint of interest in image I2.
% structA: images struct.
%   - I2x: (2 x 1 struct of R^2 matrix) gradient of image I2.
% structB: graph struct.
%   - V1: (cell of (d x 1) vectors), vertex set in image I1.
%   - V2: (cell of (d x 1) vectors), initial vertex set in image I2.
%   - E1: (cell of spline structs), edge set in image I1.
%   - E2: (cell of spline structs), initial edge set in image I2.
%   - E2new: (cell of structs) the splines in image I2, some of which are
%   adjacent to vertex v in image I2 and moved in the CURRENT ITERATION.
%   If E2new is set equal to E2, then that assumes all vertices
%   are updated simultaneously in Newton's method.
%   - adjList
% structC: parameters struct.
%   - alpha: ([0,1] scalar) percent of weight on vertex cost.
%   - w: (odd scalar) width of the vertex window in # of pixels.
% structD: functions struct.
%   - fvert
%   - fedge
%   - DfedgeDv
%   - DfvertDv
%
% OUTPUTS
% g: gradient of graph cost function wrt vertex v.
%
% CALLEE functions
%   vWeight
%   eWeight
%   vNeighborE
%
% @author Roger Zou
% @date 7/29/15

% get input variables
I2x = structA.I2x;
V1 = structB.V1;
V2 = structB.V2;
E1 = structB.E1;
E2 = structB.E2;
E2new = structB.E2new;
adjList = structB.adjList;
alpha = structC.alpha;
w = structC.w;
fvert = structD.fvert;
fedge = structD.fedge;
DfvertDv = structD.DfvertDv;
DfedgeDv = structD.DfedgeDv;
d = 2;

% get corresponding vertices in I1 and I2
u = V1{uInd};
v = V2{vInd};

Grad_v = zeros(d, 1);
Grad_e = zeros(d, 1);
g = zeros(d, 1);

if sum(isnan(u) + isnan(v))
    return
end

%% Vertex component

% (w x w) | compute weights for vertices
vweight = vWeight(structC);

% (w x w) | vertex residual at every point in window
f_k = fvert(u, v, structA, structC);

% (d x w x w) | gradient of vertex residual gradient at every pt in window (col;row)
dfvertdv = DfvertDv(v, I2x, w);

% (d x 1) | Compute vertex component of gradient
Grad_v_mat = NaN(d, w, w);
for i=1:w
    for j=1:w
        Grad_v_mat(:,i,j) = f_k(i,j) * dfvertdv(:,i,j) * vweight(i,j);
    end
end
Grad_v_mat = sum(Grad_v_mat, 3);
Grad_v_mat = sum(Grad_v_mat, 2);

% add value to vertex component
Grad_v = Grad_v + Grad_v_mat;

%% Edge component
    
% check if the edge component should be considered
if ~isempty(E1) && ~isempty(E2) && length(E1)==length(E2)

    % (1 x l) | get weight (1D gaussian kernel)
    eweights = eWeight(structC);

    % get vertices and edges neighboring u and v in images I1 and I2
    [~, edgeJ, vertJ] = vNeighborE(vInd, V2, E2, adjList);
    [~, edgeJnew, vertJnew] = vNeighborE(vInd, V2, E2new, adjList);
    [~, edgeI] = vNeighborE(uInd, V1, E1, adjList);
    if length(edgeI) ~= length(edgeJ) || length(edgeJ) ~= length(edgeJnew)
        error('DgcostDv: number of edges not equal!');
    else
        % M is the number of matching edges
        M = length(edgeJ);
    end

    % N entry cell array of (n x l) matrices | edge residuals for each relevant
    % spline
    f_kl = cell(M, 3);
    % cell array of (d x n x l) matrices | gradient of edge residuals for each
    % relevant spline. Divided into gradient of f wrt vertex v_k, the vertex of 
    % interest, and gradient wrt vertex v_l, the corresponding vertex
    % connected to vertex v_k.
    gradKf_kl = cell(M, 1);
    for ii=1:M
        
        edgei = edgeI{ii};
        edgej = edgeJ{ii};
        
        [f_ij, dsI, dxI] = fedge(edgei, edgej, structA, structC);
        f_kl(ii,:) = {f_ij, dsI, dxI};
        
        [dfedgedv1, dfedgedv2] = DfedgeDv(edgei, edgej, structA, structC, structD);
        if v == edgej.control(:,1)
            gradKf_kl{ii} = dfedgedv1;
        elseif v == edgej.control(:,end);
            gradKf_kl{ii} = dfedgedv2;
        end
    end

    % cell array of d x 1 vectors | compute v_j - v_j^{-} for each j
    diffv = cell(M, 1);
    for ii=1:M
        diffv{ii} = vertJnew{ii} - vertJ{ii};
    end

    % cell array of (d x 1) vectors of edge gradients
    edgeGrad = cell(M, 1);
    for ii=1:M
        % get (d x n x l) matrices of the edge gradient at each dx and ds
        eachf_kl = f_kl{ii,1};
        [n, l] = size(eachf_kl);
        eachgradJ = NaN(d, n, l);
        for i=1:n
            for j=1:l
                f_kl_ij = eachf_kl(i,j);
                dsI = f_kl{ii,2};
                dxI = f_kl{ii,3};
                gradKf_kl_ij = gradKf_kl{ii}(:,i,j);
                diffv_ij = diffv{ii};
                eachgradJ(:,i,j) = (sum(gradKf_kl_ij .* diffv_ij) + f_kl_ij) .* gradKf_kl_ij .* eweights(l);
            end
        end
        % sum dx and ds
        gradJ = NaN(d, 1);
        for i=1:d
            gradJ(i) = integrate(squeeze(eachgradJ(i,:,:)), dsI, dxI);
        end
        edgeGrad{ii} = gradJ;
    end
    
    % sum over each edge gradient component for overall edge gradient
    for ii=1:M
        Grad_e = Grad_e + edgeGrad{ii};
    end
    
end

%% Sum over all components

g = alpha * Grad_v + (1-alpha) * Grad_e;

end
