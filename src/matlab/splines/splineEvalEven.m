% Evaluate the curve and possibly spline matrices, tangents, normals, and
% needles for spline s, and ensure an approximately uniform sampling along
% the curve. See splineEval for more information.

function s = splineEvalEven(s, makeP, makeT, makeNeedles)

if nargin < 2 || isempty(makeP)
    makeP = true;
end

if nargin < 3 || isempty(makeT)
    makeT = true;
end

if nargin < 4 || isempty(makeNeedles)
    makeNeedles = true;
end

s = splineEval(s, false, false, false);
for i = 1:3
    s.t = reparm(s.t);
    s = splineEval(s, false, false, false);
end
s.t = reparm(s.t);
s = splineEval(s, makeP, makeT, makeNeedles);

    % Reparameterize to arclength (approximately)
    function y = reparm(x)
        c = s.curve;
        a = [0 cumsum(sqrt(diff(c(1, :)).^2 + diff(c(2, :)).^2))];
        t = linspace(0, a(end), length(s.t));
        y = interp1(a, x, t);
    end

end