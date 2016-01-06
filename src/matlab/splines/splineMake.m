% Make a spline with the given control points (a dim by n matrix), order,
% and multiplicities. If open is true, the spline is open; otherwise, it is
% periodic. The parameter n is the number of points to draw (can be 0).
%
% @author Carlo Tomasi and Roger Zou
% @date 5/20/15

function s = splineMake(control, order, mult, open, n, makeNeedles)

if nargin < 6 || isempty(makeNeedles)
    makeNeedles = true;
end

% save initialization parameters
s.control = control;
s.order = order;
s.mult = mult;
s.open = open;  % Remember if the spline is open (true) or periodic (open false)
s.n = n;
s.d = size(control, 1);
s.k = size(control, 2)-2;

[s.knot, s.Ncv, s.kIdxRange] = knots(order, open, mult);

if size(s.control, 2) ~= s.Ncv
    error('The number of control points is inconsistent with the given knot multiplicities')
end

s.matrix = CoxDeBoor(order, s.knot, open);

% Curve parameter values
s.t = linspace(s.knot(s.kIdxRange(1)), s.knot(s.kIdxRange(2)), n);

if isempty(s.t)
    s.curve = [];
else
    if ~open
       % Add in the knot values in the drawable range
        s.t = unique([s.t s.knot(s.kIdxRange(1)):s.knot(s.kIdxRange(2))]);
    end
    
    % Compute points on the spline, spline matrices, tangents, normals, and needles
    s = splineEvalEven(s, true, true, makeNeedles);
end