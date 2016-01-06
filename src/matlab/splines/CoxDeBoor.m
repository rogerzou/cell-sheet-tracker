% Make the B-spline basis matrices for the given order and knot vector
% by the method introduced in
%    Kaihuai Qin, "General matrix representations for B-splines,"
%    The Visual Computer, 16:177-186, Springer, 2000.
%
% The spline order defaults to 4. If knot is unspecified or empty, the
% matrix for the Bezier spline of the given order results.
%
% The function also returns the knot vector. This is useful if the vector
% is generated internally rather than passed as an input argument.
% Use the function knots.m to generate open or periodic knot vectors.

function [M, knot] = CoxDeBoor(order, knot, open)

if nargin < 1 || isempty(order)
    order = 4;    % Cubic spline
elseif order < 1
    error('Order must be at least 1');
end

if nargin < 2 || isempty(knot)
    knot = knots(order);    % Knot vector for the Bezier spline
end

% Find indices of span starts
span = spans(knot, order, open);

% Number of spans
L = length(span);

if L < 1
    error('Need at least one span')
end

% Top recursive call of Cox-DeBoor
M = cdb(order);

    % Recursive Cox-DeBoor routine
    function B = cdb(k)
        if k == 1
            B = ones(1, 1, L);
        else
            C = cdb(k - 1);
            z = zeros(k - 1, 1);
            B = zeros(k, k, L);
            for s = 1:L
                i = span(s);
                j = (i - k + 2):i;
                num0 = knot(i) - knot(j);
                num1 = (knot(i+1) - knot(i)) * ones(1, length(j));
                denom = knot(j + k - 1) - knot(j);
                nzd = denom ~= 0;
                denom = denom(nzd);
                d0 = zeros(1, length(j));
                d0(nzd) = num0(nzd) ./ denom;
                d1 = zeros(1, length(j));
                d1(nzd) = num1(nzd) ./ denom;
                U = [diag(1 - d0) z] + [z diag(d0)];
                V = [diag(-d1) z] + [z diag(d1)];
                B(:, :, s) = [C(:, :, s); z'] * U + [z'; C(:, :, s)] * V;
            end
        end
    end
end
