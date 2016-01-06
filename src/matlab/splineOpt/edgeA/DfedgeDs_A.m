function [ DfedgeDs ] = DfedgeDs_A( spline, Jx, l )
%DFEDGEDS_A DfedgeDq for free splines.
%
% INPUTS
% spline: (struct) spline.
% Jx: (2x1 cell) gradient of image at every point.
% l: (odd scalar) spread of the spline edge in # of pixels.
%
% OUTPUTS
% DfedgeDs
%
% CALLEE functions
%   nablaJ
%   DyDqT_A
%
% @author Roger Zou
% @date 6/9/15

% get gradient of J evaluated at spline points (n x l x d)
nablaj = nablaJ(spline, Jx, l);

% get component of the gradient of f at each sample point (d x r x n x l)
DyDqt = DyDqT_A(spline, l);
[~, r, n, l] = size(DyDqt);

% compute derivative
nablaj = permute(nablaj, [3,1,2]);
DfedgeDs = NaN(n, l, r);
for i=1:n
    for j=1:l
        dydq_ij = DyDqt(:,:,i,j)';
        nablaj_ij = nablaj(:,i,j);
        DfedgeDs(i,j,:) = - dydq_ij * nablaj_ij;
    end
end

end
