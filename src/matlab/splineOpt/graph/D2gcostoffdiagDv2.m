function [ Gr ] = D2gcostoffdiagDv2( u1Ind, u2Ind, v1Ind, v2Ind, structA, structB, structC, structD )
%D2GCOSTOFFDIAGDV2 An off-diagional component of the graph cost Gramian wrt
%vertices v_1 and v_2.
%
% INPUTS
% u1Ind: (integer) index of endpoint 1 of interest in image I1.
% u2Ind: (integer) index of endpoint 2 of interest in image I1.
% v1Ind: (integer) index of endpoint 1 of interest in image I2.
% v2Ind: (integer) index of endpoint 2 of interest in image I2.
% structA: images struct.
%   - I1: (R^2 matrix) image I1, the previous image in a sequence.
% structB: graph struct.
%   - V1: (cell of (d x 1) vectors), vertex set in image I1.
%   - V2: (cell of (d x 1) vectors), initial vertex set in image I2.
%   - E1: (cell of spline structs), edge set in image I1.
%   - E2: (cell of spline structs), initial edge set in image I2.
%   - adjList
% structC: parameters struct.
%   - alpha: ([0,1] scalar) percent of weight on vertex cost.
% structD: functions struct.
%   - DfedgeDv
%
% OUTPUTS
% Gr: off-diagonal component of the graph cost function Gramian wrt
% vertices k and l
%
% CALLEE functions
%   eWeight
%   splineValue
% 
% @author Roger Zou
% @date 7/29/15

% get input variables
I1 = structA.I1;
V1 = structB.V1;
V2 = structB.V2;
E1 = structB.E1;
E2 = structB.E2;
adjList = structB.adjList;
alpha = structC.alpha;
DfedgeDv = structD.DfedgeDv;
d = 2;

u1 = V1{u1Ind};
v1 = V2{v1Ind};
u2 = V1{u2Ind};
v2 = V2{v2Ind};

Gram_e = zeros(d);
Gr = zeros(d);

if sum(isnan(u1) + isnan(u2) + isnan(v1) + isnan(v2))
    return
end

%% Edge component of Gramian

% check if the edge component should be considered
if ~isempty(E1) && ~isempty(E2) && length(E1)==length(E2)

    % (1 x l) | get weight (1D gaussian kernel)
    eweights = eWeight(structC);

    % (d x n x l) matrix | gradient of edge residuals for spline of interest.
    EdgesU = uvEdge(u1Ind, u2Ind, V1, E1, adjList);
    EdgesV = uvEdge(v1Ind, v2Ind, V1, E2, adjList);
    if isempty(EdgesU) && isempty(EdgesV)
        return;
    end

    % iterate over each edge between two vertices
    for ii=1:numel(EdgesV)

        % get edge of interest in both image I1 and I2
        edgeu = EdgesU{ii};
        edgev = EdgesV{ii};

        % compute edge residual
        [dfedgedv1, dfedgedv2] = DfedgeDv(edgeu, edgev, ...
            structA, structC, structD);
        if v1 == edgev.control(:,1)
            gradKf_kl = dfedgedv1;
            gradLf_kl = dfedgedv2;
        elseif v1 == edgev.control(:,end)
            gradKf_kl = dfedgedv2;
            gradLf_kl = dfedgedv1;
        end

        % get sizes
        [d, n, l] = size(gradKf_kl);

        % get ds and dx
        [~, dsI, dxI] = splineValue(edgeu, I1, l);

        % get (d x d x n x l) matrix of the edge Gramian at each dx and ds
        eachgramJ = NaN(d, d, n, l);
        for jj=1:n
            for kk=1:l
                gradKf_kl_ij = gradKf_kl(:,jj,kk);
                gradLf_kl_ij = gradLf_kl(:,jj,kk);
                eachgramJ(:,:,jj,kk) = gradKf_kl_ij * gradLf_kl_ij' .* eweights(l);
            end
        end
        % sum dx and ds
        Gram_e = NaN(d, d);
        for jj=1:d
            for kk=1:d
                Gram_e(jj,kk) = integrate(squeeze(eachgramJ(jj,kk,:,:)), dsI, dxI);
            end
        end
    end

end

%% Sum over all components (just edge)

Gr = (1-alpha) * Gram_e;

end
