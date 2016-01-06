function [ g, Gr ] = D2ecostDs2( sU, sV, structA, structC, structD )
%D2ECOSTDS2 Derivative of edge cost with respect to independent variables.
% NOTE: nWidth MUST be odd, and the gaussian weight supports 2.5*SD.
%
% INPUTS
% sU: (struct) spline in image I.
% sV: (struct) spline in image J.
% structA: images struct.
%   - I2x: (2 x 1 cell of R^2 matrix) gradient of image I2.
% structC: parameters struct.
%   - l: (odd scalar) spread of the spline edge in # of pixels.
% structD: functions struct.
%   - fedge: residual at every point along spline.
%   - DfedgeDs: gradient of residual at every point along spline.
%
% OUTPUTS
% g: gradient of edge cost function wrt spline independent variables.
% Gr: Gramian of edge cost function wrt spline independent variables.
%
% CALLEE functions
%   eWeight
%   integrate
%
% @author Roger Zou
% @date 7/29/15

% get input variables
I2x = structA.I2x;
l = structC.l;
fedge = structD.fedge;
DfedgeDs = structD.DfedgeDs;

% validate inputs
if ~mod(l, 2)
    error('D2ECOSTDS2: nWidth must be odd!');
end

% weight (1D gaussian kernel)
eweight = eWeight(structC);

% edge residual
[f_e, dsI, dxI] = fedge(sU, sV, structA, structC);

% gradient of edge residual
df_e = DfedgeDs(sV, I2x, l);
[n, l, r] = size(df_e);

% compute gradient and Gramian at each point
df_e = permute(df_e, [3,1,2]);
g_ce = NaN(n, l, r);
Gr_ce = NaN(n, l, r, r);
for i=1:n
    for j=1:l
        f_ij = f_e(i,j);
        df_ij = df_e(:,i,j);
        Gr_ce(i,j,:,:) = df_ij * df_ij';
        g_ce(i,j,:) = f_ij * df_ij;
    end
end

% sum to obtain gradient and Gramian
g = NaN(r, 1);
Gr = NaN(r, r);
for i=1:r
    g_i = g_ce(:,:,i);
    g_i = bsxfun(@times, g_i, eweight);
    g(i) = integrate(g_i, dsI, dxI);
    for j=1:r
        Gr_ij = Gr_ce(:,:,i,j);
        Gr_ij = bsxfun(@times, Gr_ij, eweight);
        Gr(i,j) = integrate(Gr_ij, dsI, dxI);
    end
end

end
