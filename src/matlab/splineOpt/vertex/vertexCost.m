function c = vertexCost(u, v, structA, structC, structD)

fvert = structD.fvert;

vweight = vWeight(structC);
f_i = fvert(u, v, structA, structC);
costV = double(f_i).^2 .* vweight;
c = sum(costV(:));
end