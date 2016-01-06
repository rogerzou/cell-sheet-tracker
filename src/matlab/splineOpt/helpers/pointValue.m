function [ Values, A] = pointValue( v, C, w )
%POINTVALUE Evaluate point v in square with length w in cost function C
%
% INPUTS
% v: (d x 1) vector of the specific point of interest in (col;row) format
% C: (R^2 matrix) image cost function.
% w: (odd scalar) side length of the square around the v to evaluate C
%
% OUTPUTS
% Values: (w x w) cost matrix, entry returns NaN if out of bounds.
% A: (w x w x d) coordinate matrix of the evaluated points in C in (col;row)
% 
% @author Roger Zou
% @date 3/20/15

% validate inputs
if ~mod(w, 2)
    error('POINTVALUE: w must be odd!');
end

C = double(C);

% compute window coordinates
halfw = floor(w/2);
[v_c, v_r] = meshgrid(-halfw:halfw);
AC = (v_c + v(1));
AR = (v_r + v(2));
A = cat(3, AC, AR);

% get sizes
[rmax, cmax] = size(C);
[rT, cT] = size(AC);

% check out of bounds
if min(AC(:)) < 1 || min(AR(:)) < 1 || max(AC(:)) > cmax || max(AR(:)) > rmax
    error('POINTVALUE: out of bounds!');
end

% compute values at points A in cost function C
Values = interp2(C, AC(:), AR(:), 'linear');
Values = reshape(Values, [rT, cT]);

end
