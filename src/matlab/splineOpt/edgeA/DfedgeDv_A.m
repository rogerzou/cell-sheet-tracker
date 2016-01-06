function [ DfedgeDa, DfedgeDb ] = DfedgeDv_A( sV1, sV2, structA, structC, structD )
%DFEDGEDV_A Compute the gradient of the edge residual function f
%
% INPUTS
% sV1: (struct) the spline in image I1. This spline is fixed, and used to
% compute the optimal sV interior control points.
% sV2: (struct) the spline in image I2. The optimality condition on the
% interior control points is a measure of similarity between sV1 in I1 and
% sV2 in I2.
% structA: images struct.
%   - I2x: (2 x 1 cell of R^2 matrix) gradient of image I2.
% structC: parameters struct.
%   - l: (odd scalar) width of spline in # of pixels.
% structD: functions struct.
%
% OUTPUTS
% DfedgeDa: (d x n x l) the gradient of the residual function f for
% edges wrt endpoint vertex 1 evaluated at each point along spline.
% DfedgeDb: (d x n x l) the gradient of the residual function f for
% edges wrt endpoint vertex 2 evaluated at each point along spline.
%
% CALLEE functions
%   DyDqT_A
%   D2ecostDs2
%   DqDvT_A
%
% @author Roger Zou
% @date 6/12/15

% get input variables
I2x = structA.I2x;
l = structC.l;
d = 2;
n = sV2.n;

% get gradient of J evaluated at spline points. nablan is (d x n x l)
nablaj = nablaJ(sV2, I2x, l);

% get dydqt, a (d x m x n x l) matrix for each point i \in 1:n, j \in 1:l
dydqt = DyDqT_A(sV2, l);

% get jacobians of xi wrt endpoints a and b, both (m x d) matrices.
[~, Gr_e] = D2ecostDs2(sV1, sV2, structA, structC, structD);
[dqdat, dqdbt] = DqDvT_A(Gr_e);

% construct derivative of f for edges wrt each endpoints a and b.
% (d x n x l)
nablaj = permute(nablaj, [3,1,2]);
DfedgeDa = NaN(d, n, l);
DfedgeDb = NaN(d, n, l);
for i=1:n
    for j=1:l
        dydqt_ij = dydqt(:,:,i,j)';
        djdy = nablaj(:,i,j);
        DfedgeDa(:,i,j) = - dqdat' * dydqt_ij * djdy;
        DfedgeDb(:,i,j) = - dqdbt' * dydqt_ij * djdy;
    end
end

end
