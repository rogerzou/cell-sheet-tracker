function [g, Gr] = vertexOpt(u, v, structA, structC, structD)

I2x = structA.I2x;
w = structC.w;
fvert = structD.fvert;
DfvertDv = structD.DfvertDv;
d = 2;

% gradient
% (w x w) | compute weights for vertices
vweight = vWeight(structC);
% (w x w) | vertex residual at every point in window
f_k = fvert(u, v, structA, structC);
% (d x w x w) | gradient of vertex residual gradient at every pt in window (col;row)
dfvertdv = DfvertDv(v, I2x, w);
% (d x 1) | Compute vertex component of gradient
Grad_v_mat = NaN(d, w, w);
for ii=1:w
    for jj=1:w
        Grad_v_mat(:,ii,jj) = f_k(ii,jj) * dfvertdv(:,ii,jj) * vweight(ii,jj);
    end
end
g = sum(sum(Grad_v_mat, 3), 2);

% Gramian
% (d x w x w) | compute vertex residual gradient of square centered at 
% vertex point v
dfvertdv = DfvertDv(v, I2x, w);

% (w x w) | compute weights for vertices
vweight = vWeight(structC);

% (d x d) | Compute vertex component of Gramian
Gram_v_mat = NaN(d, d, w, w);
for ii=1:w
    for jj=1:w
        Gram_v_mat(:,:,ii,jj) = dfvertdv(:,ii,jj) * dfvertdv(:,ii,jj)' * vweight(ii,jj);
    end
end
Gr = sum(sum(Gram_v_mat, 4), 3);

end