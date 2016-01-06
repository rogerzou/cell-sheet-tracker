% Compute a dim by m matrix spl.curve whose columns are the values of the given
% B-spline spl (made with splineMake) for the parameter values in spl.t .
% The columns of spl.tangent are unit vectors along the tangent to the spline
% curve.
% If makeP is true (it is by default), also compute a (dim * m) by
% (dim * n) matrix spl.P such that if c = spl.curve(:) then c = spl.P * v where
% v = spl.control(:) . Here, spl.control is a dim by n matrix of B-spline
% control points.
% If makeT is true (it is by default), also compute a (dim * m) by
% (dim * n) matrix spl.T such that if Tg = spl.tangent(:) then Tg = spl.T * v where
% v = spl.control(:) .
% If makeT is true and dim == 2 (spline on the plane), also compute unit
% normal vectors spl.normal and a (dim * m) by
% (dim * n) matrix spl.N such that if No = spl.normal(:) then No = spl.N * v where
% v = spl.control(:) .
% If makeNeedles, also computes the needles and adds spl.needle

function spl = splineEval(spl, makeP, makeT, makeNeedles)

if nargin < 2 || isempty(makeP)
    makeP = true;
end

if nargin < 3 || isempty(makeT)
    makeT = true;
end

if nargin < 4 || isempty(makeNeedles)
    makeNeedles = true;
end

M = spl.matrix;
order = size(M, 1);
knot = spl.knot;
V = spl.control;
t = spl.t;
open = spl.open;
tmin = knot(spl.kIdxRange(1));
tmax = knot(spl.kIdxRange(2));
trange =  tmax - tmin;

% Check that we have the correct number of points
if open
    Ncv = length(knot) - order;
else
    Ncv = length(knot) - order - 1;
end
if size(V, 2) ~= Ncv
    error('The number of control points is inconsistent with the spline matrices')
end
dim = size(V, 1);
V = V';

m = length(t);

% Find indices of span starts
span = spans(knot, order, open);

L = length(span);
if L ~= size(M, 3)
    error('The number of spline matrices is different from the number of spans')
end

if makeP
    n = size(spl.control, 2);
    P = zeros(m, n); % We'll make it bigger later with a Kronecker product
end

if makeT
    n = size(spl.control, 2);
    T = zeros(m, n); % We'll make it bigger later with a Kronecker product
end

% Evaluate the spline
C = zeros(m, dim);
Tg = zeros(m, dim);
for j = 1:m
    if open
        tj = t(j);
    else
        tj = mod(t(j) - tmin, trange) + tmin;
    end
    
    % Which span tj is in
    si = find(knot(span) <= tj & knot(span + 1) > tj);
    if open
        if isempty(si)
            si = length(span);
        end
    end
    s = span(si);
        
    % Knot value at the span endpoints
    ti = knot(s);
    tip = knot(s + 1);
    
    % Normalized span coordinate
    u = (tj - ti) / (tip - ti);
    
    % Coordinate power vector
    p = powers(u, order);
    
    % Extract control points for this span
    rs = mmod(s - order + 1, Ncv);
    es = mmod(rs:(rs + order - 1), Ncv);
    Vspan = V(es, :);
    
    % Compute the point on the spline
    pM = p * M(:, :, mmod(si, L));
    C(j, :) = pM * Vspan;
    
    % Compute the unit tangent vectors
    dp = dpowers(u, order);
    tM = dp * M(:, :, mmod(si, L));
    Tg(j, :) = tM * Vspan;
    Tnorm = norm(Tg(j, :));
    Tg(j, :) = Tg(j, :) / Tnorm;
    tM = tM / Tnorm;

    if makeP
        E = zeros(order, n);
        E(:, es) = eye(order);
        P(j, :) = pM * E;
    end
    
    if makeT
        E = zeros(order, n);
        E(:, es) = eye(order);
        T(j, :) = tM * E;
    end
end

spl.curve = C';
spl.tangent = Tg';

if makeP
    spl.P = kron(P, eye(dim));
end

if makeT
    spl.T = kron(T, eye(dim));
    if dim == 2 % Also make a matrix to compute normals from control points
        spl.N = zeros(size(spl.T));
        spl.N(1:2:end, :) = -spl.T(2:2:end, :);
        spl.N(2:2:end, :) = spl.T(1:2:end, :);
        spl.normal = reshape(spl.N * spl.control(:), size(Tg'));
    end
end

if makeNeedles
    spl = splineNeedles(spl);
end

% First kk powers of the scalar x (starting at x^0)
    function q = powers(x, kk)
        q = zeros(1, kk);
        q(1) = 1;
        for jj = 2:kk
            q(jj) = x * q(jj-1);
        end
    end

% Derivatives of the first kk powers of the scalar x
    function q = dpowers(x, kk)
        q = zeros(1, kk);
        q(2) = 1;
        for jj = 3:kk
            q(jj) = (jj - 1) / (jj - 2) * x * q(jj-1);
        end
    end
end