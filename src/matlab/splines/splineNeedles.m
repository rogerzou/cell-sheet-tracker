% Add needles to spline s, extending up to d units in either direction.
%
% A needle at point p on s.curve is a line segment orthogonal to s.curve at
% p. It extends d units outwards and min(d, w) units inwards, where w is
% the distance between p and the skeleton of the spline curve.
%
% If d is unspecified, it is set to the maximum value of d over the points
% that define s.curve.

function s = splineNeedles(s, d)

if nargin < 3
    d = [];
end

p = s.curve';
n = size(p, 1);

% Remember state of Delaunay warning, and turn it off
msgid = 'MATLAB:DelaunayTri:DupPtsConsUpdatedWarnId';
wst = warning('query', msgid);
warning('off', msgid);

% Make a triangulation of the interior of the spline, and find triangle
% circumcenters. These are on the skeleton.
constraint = [(1:(n-1))' (2:n)'; n 1];
dt = DelaunayTri(p(:, 1), p(:, 2), constraint);
inside = dt.inOutStatus();
tr = TriRep(dt(inside, :), dt.X);
cc = tr.circumcenters();

% Distances between spline curve points and nearest circumcenters
dtc = DelaunayTri(cc(:, 1), cc(:, 2));
pi = nearestNeighbor(dtc, p);
delta = p - cc(pi, :);
w = sqrt(sum(delta.^2, 2));

% Restore original warning state
warning(wst.state, msgid);

% Needle unit vectors
s.normal = s.tangent;
s.normal = [s.normal(2, :); -s.normal(1, :)];

% Orient normal vectors so they point outwards
small = norm(min(p(2:end, :) - p(1:(end-1), :)))/10;
q = mean(p(1:2, :) + small * s.normal(:, 1:2)', 1);
if inpolygon(q(1), q(2), p(:, 1), p(:, 2))
    s.normal = -s.normal;
end

if isempty(d)
    d = max(w);
end

% Add the needles
p = p';
s.needle.in = p - s.normal .* ([1; 1] * w');
s.needle.out = p + s.normal * d;

