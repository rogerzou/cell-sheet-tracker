function [ g, Gr ] = D2gcostDv2( structA, structB, structC, structD, parallel )
%D2GCOSTDV2 Gradient and Gramian of the graph cost function.
%
% INPUTS
% structA: images struct.
% structB: graph struct.
%   - V2: (cell of (d x 1) vectors), vertex set in image I2.
% structC: parameters struct.
% structD: functions struct.
% parallel: (boolean, default false) parallelization toggle.
%
% OUTPUTS
% g: gradient of graph cost function wrt vertex set V.
% Gr: Gramian of the graph cost function wrt vertex set V.
%
% CALLEE functions
%   DgcostDv
%   D2gcostdiagDv2
%   D2gcostoffdiagDv2
%
% @author Roger Zou
% @date 7/29/15

% check for parallelization
if nargin < nargin('D2gcostDv2') || isempty(parallel)
    parallel = false;
end

% get input variables
V2 = structB.V2;
N = length(V2);
d = 2;
n = N*d;

%% Gradient (PARFOR POSSIBLE)
if parallel
    g = zeros(d, N);
    parfor ii=1:N
        g(:,ii) = DgcostDv(ii, ii, structA, structB, structC, structD);
    end
    g = g(:);
else
    g = zeros(n, 1);
    r1 = 1 + d * ( (1:N) - 1 );
    r2 = d * (1:N);
    for ii=1:N
        g(r1(ii):r2(ii)) = DgcostDv(ii, ii, structA, structB, structC, structD);
    end
end

%% Gramian (PARFOR POSSIBLE)
if parallel
    triuInd = triu(true(N,N));
    triuInd = triuInd(:);
    x = (1:length(triuInd))';
    [R, C] = ind2sub([N,N], x(triuInd));
    M = length(R);
    GrI = cell(M,1);
    parfor ii=1:M
        ri = R(ii);
        ci = C(ii);
        if ri == ci     % diagonal component of Gramian
            Gr_ij = D2gcostdiagDv2( ri, ri, ...
                structA, structB, structC, structD);
        else
            Gr_ij = D2gcostoffdiagDv2( ri, ci, ri, ci, ...
                structA, structB, structC, structD);
        end
        GrI{ii} = Gr_ij;
    end
    Gr = zeros(n, n);
	r1 = 1 + d * ( R - 1 );
    r2 = d * R;
    c1 = 1 + d * ( C - 1 );
    c2 = d * C;
    for ii=1:M
        Gr(r1(ii):r2(ii), c1(ii):c2(ii)) = GrI{ii};
    end
    Gr = triu(Gr) + triu(Gr, 1)';
else
    Gr = zeros(n, n);
    r1 = 1 + d * ( (1:N) - 1 );
    r2 = d * (1:N);
    c1 = 1 + d * ( (1:N) - 1 );
    c2 = d * (1:N);
    for ii=1:N
        for jj=1:N
            if ii == jj     % diagonal component of Gramian
                Gr_ij = D2gcostdiagDv2( ii, ii, ...
                    structA, structB, structC, structD);
            else            % off-diagonal component of Gramian
                Gr_ij = D2gcostoffdiagDv2( ii, jj, ii, jj, ...
                    structA, structB, structC, structD);
            end
            % insert Gramian component into Gramian
            Gr(r1(ii):r2(ii), c1(jj):c2(jj)) = Gr_ij;
        end
    end
end

