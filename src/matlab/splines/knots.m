% Make a B-spline knot vector of the given order and with internal
% breakpoints having the multiplicities given in the vector intmult.
% If open is true (the default) the spline is open, and is periodic
% otherwise.
%
% Multiple breakpoints, if any, are the only departure from uniformity.
% The default order is 4 (cubic spline), and there are no internal
% multiplicities by default.
%
% When open is true and inmult is empty, i.e.,
%    knot = knots(order);
% the knot vector for a Bezier spline of the given order is returned.
% When open is false and inmult is empty, i.e.,
%    knot = knots(order, false);
% the knot vector for a uniform periodic spline of the given order is returned.
%
% The function also returns the required number Ncv of control vertices, the
% smallest and largest allowable parameter value (tRange) and the indices
% of the corresponding knot values (kIdxRange).
% Control vertices need not be distinct.

function [knot, Ncv, kIdxRange] = knots(order, open, intmult)

if nargin < 1 || isempty(order)
    order = 4;    % Cubic spline
elseif order < 1
    error('Order must be at least 1')
end

if nargin < 2 || isempty(open)
    open = 'true';
end

if nargin < 3 || isempty(intmult)
    intmult = [];    % A single span with no internal breakpoints
% This gives the Bezier spline if open is true, and the uniform periodic
% spline otherwise
end

if open
    % Assign appropriate endpoint mutliplicities
    mult = [order intmult order];
else % Periodic spline
    o = ones(1, order);
    mult = [o intmult o];
end

Nbrk = length(mult);    % Number of breakpoints
L = Nbrk - 1;           % Number of spans
M = sum(mult - 1);      % Multiple knot count
Nk = L + M + 1;         % Number of knots

% Make knots with the given multiplicities
knot = zeros(1, Nk);
rs = 1;
for bp = 1:Nbrk
    re = rs + mult(bp) - 1;
    knot(rs:re) = bp - 1;
    rs = re + 1;
end

if ~open
    knot = knot - order + 1;
end

if open
    Ncv = Nk - order;
else
    Ncv = Nk - order - 1;
end

% Range of allowable curve parameter values and correspoinding knot
% indices. Values outside this range are for clamping the endpoints or
% closing the curve
if open
    kIdxRange = [order length(knot) - order + 1];
else
    % Wider parameter range to close the curve
    kIdxRange = [order length(knot) - order + 2];
end

