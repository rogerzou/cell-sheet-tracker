function [ c, cV, cE ] = gcost(structA, structB, structC, structD, parallel)
%GCOST Graph cost.
%
% INPUTS
% structA: images struct.
% structB: graph struct.
%   - V1: (cell of (d x 1) vectors), vertex set in image I1.
%   - V2: (cell of (d x 1) vectors), initial vertex set in image I2.
%   - E1: (cell of spline structs), edge set in image I1.
%   - E2: (cell of spline structs), initial edge set in image I2.
% structC: parameters struct.
%   - alpha: weight of vertex cost (1 = all vertex, 0 = all edge)
% structD: functions struct.
%   - fvert
%   - fedge
% parallel: (boolean, default false) parallelization toggle.
%
% OUTPUTS
% c: graph cost.
% cV: unweighted vertex cost component.
% cE: unweighted edge cost component.
%
% CALLEE functions
%   vWeight
%   eWeight
%   integrate
%
% @author Roger Zou
% @date 7/29/15

% get input parameters
V1 = structB.V1;
V2 = structB.V2;
E1 = structB.E1;
E2 = structB.E2;
alpha = structC.alpha;
fvert = structD.fvert;
fedge = structD.fedge;
N = numel(V2);
M = numel(E2);

% vertex weights (gaussian)
vweight = vWeight(structC);

% edge weights (gaussian)
eweight = eWeight(structC);


%% compute vertex cost
cV = zeros(N,1);
if parallel
    parfor ii=1:N
        u = V1{ii};
        v = V2{ii};
        if ~sum(isnan(u) + isnan(v))
            f_i = fvert(u, v, structA, structC);
            costV = double(f_i).^2 .* vweight;
            cV(ii) = sum(costV(:));
        end
    end
else
    for ii=1:N
        u = V1{ii};
        v = V2{ii};
        if ~sum(isnan(u) + isnan(v))
            f_i = fvert(u, v, structA, structC);
            costV = double(f_i).^2 .* vweight;
            cV(ii) = sum(costV(:));
        end
    end
end
cV = sum(cV);

%% compute edge cost
cE = zeros(M,1);
if parallel
    parfor ii=1:M
        % get splines and cost
        uspline = E1{ii};
        vspline = E2{ii};
        if ~isempty(uspline) && ~isempty(vspline)
            [f_ij, dsI, dxI] = fedge(uspline, vspline, structA, structC);
            costE = bsxfun(@times, (f_ij).^2, eweight);
            cE(ii) = integrate(costE, dsI, dxI);
        end
    end
else
    for ii=1:M
        % get splines and cost
        uspline = E1{ii};
        vspline = E2{ii};
        if ~isempty(uspline) && ~isempty(vspline)
            [f_ij, dsI, dxI] = fedge(uspline, vspline, structA, structC);
            costE = bsxfun(@times, (f_ij).^2, eweight);
            cE(ii) = integrate(costE, dsI, dxI);
        end
    end
end
cE = sum(cE);

%% compute graph cost
c = alpha * cV + (1-alpha) * cE;
