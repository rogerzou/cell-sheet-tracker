function [ Gr ] = D2gcostdiagDv2( uInd, vInd, structA, structB, structC, structD )
%D2GCOSTOFFDIAGDV2 A diagional component of the graph cost Gramian wrt
%vertex v.
%
% INPUTS
% uInd: (integer) index of the endpoint of interest in image I.
% vInd: (integer) index of the endpoint of interest in image J.
% structA: images struct.
%   - I1: (R^2 matrix) image I1, the previous image in a sequence.
%   - I2x: (2 x 1 struct of R^2 matrix) gradient of image I2.
% structB: graph struct.
%   - V1: (cell of (d x 1) vectors), vertex set in image I1.
%   - V2: (cell of (d x 1) vectors), initial vertex set in image I2.
%   - E1: (cell of spline structs), edge set in image I1.
%   - E2: (cell of spline structs), initial edge set in image I2.
%   - adjList
% structC: parameters struct.
%   - alpha: ([0,1] scalar) percent of weight on vertex cost.
%   - w: (odd scalar) width of the vertex window in # of pixels.
% structD: functions struct.
%   - DfvertDv
%   - DfedgeDv
%
% OUTPUTS
% Gr: diagonal component of the graph cost function Gramian wrt vertex v.
%
% CALLEE functions
%   vWeight
%   eWeight
%   vNeighborE
%   splineValue
%   integrate
%
% @author Roger Zou
% @date 7/29/15

% get input variables
I1 = structA.I1;
I2x = structA.I2x;
V1 = structB.V1;
V2 = structB.V2;
E1 = structB.E1;
E2 = structB.E2;
adjList = structB.adjList;
alpha = structC.alpha;
w = structC.w;
DfvertDv = structD.DfvertDv;
DfedgeDv = structD.DfedgeDv;
d = 2;

% get vertex in J
v = V2{vInd};

Gram_v = zeros(d);
Gram_e = zeros(d);
Gr = zeros(d);

if sum(isnan(v))
    return
end

%% Vertex component

% (d x w x w) | compute vertex residual gradient of square centered at 
% vertex point v
dfvertdv = DfvertDv(v, I2x, w);

% (w x w) | compute weights for vertices
vweight = vWeight(structC);

% (d x d) | Compute vertex component of Gramian
Gram_v_mat = NaN(d, d, w, w);
for i=1:w
    for j=1:w
        Gram_v_mat(:,:,i,j) = dfvertdv(:,i,j) * dfvertdv(:,i,j)' * vweight(i,j);
    end
end
Gram_v_mat = sum(Gram_v_mat, 4);
Gram_v_mat = sum(Gram_v_mat, 3);

Gram_v = Gram_v + Gram_v_mat;

%% Edge component

% check if the edge component should be considered
if ~isempty(E1) && ~isempty(E2) && length(E1)==length(E2)

    % (1 x l) | get weight (1D gaussian kernel)
    eweights = eWeight(structC);

    % get vertices and edges neighboring u and v in images I and J
    [~, edgeJ] = vNeighborE(vInd, V2, E2, adjList);
    [~, edgeI] = vNeighborE(uInd, V1, E1, adjList);
    if length(edgeI) ~= length(edgeJ)
        error('D2GCOSTDIAGDV2: number of edges not equal!');
    else
        % M is the number of matching edges
        M = length(edgeJ);
    end

    % cell array of (d x n x l) matrices | gradient of edge residuals for each
    % relevant spline. Divided into gradient of f wrt vertex v_k, the vertex of 
    % interest, and gradient wrt vertex v_l, the corresponding vertex
    % connected to vertex v_k.
    gradKf_kl = cell(M, 1);
    for ii=1:M
        edgei = edgeI{ii};
        edgej = edgeJ{ii};
        [dfedgedv1, dfedgedv2] = DfedgeDv(edgei, edgej, structA, structC, structD);
        if v == edgej.control(:,1)
            gradKf_kl{ii} = dfedgedv1;
        elseif v == edgej.control(:,end)
            gradKf_kl{ii} = dfedgedv2;
        end
    end

    % cell array of (d x d) matrix of edge Gramians
    edgeGram = cell(M, 1);
    for ii=1:M
        % get sizes (assuming gradKf_kl is not empty)
        [d, n, l] = size(gradKf_kl{ii});
        % get (d x d x n x l) matrices of the edge Gramian at each dx and ds
        eachgramJ = NaN(d, d, n, l);
        for i=1:n
            for j=1:l
                gradKf_kl_ij = gradKf_kl{ii}(:,i,j);
                eachgramJ(:,:,i,j) = gradKf_kl_ij * gradKf_kl_ij' .* eweights(l);
            end
        end
        % sum dx and ds
        gramJ = NaN(d, d);
        [~, dsI, dxI] = splineValue(edgeI{ii}, I1, l);
        for i=1:d
            for j=1:d
                gramJ(i,j) = integrate(squeeze(eachgramJ(i,j,:,:)), dsI, dxI);
            end
        end
        edgeGram{ii} = gramJ;
    end

    % sum over each edge gradient component for overall edge gradient
    for ii=1:M
        Gram_e = Gram_e + edgeGram{ii};
    end

end

%% Sum over all components

Gr = alpha * Gram_v + (1-alpha) * Gram_e;

end
