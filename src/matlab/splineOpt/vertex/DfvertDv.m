function [ DfvertDv ] = DfvertDv( v, Jx, w )
%DFVERTDV Compute the gradient of the vertex residual function f
%
% INPUTS
% v: (d x 1) vector of the specific point of interest in (col;row)
% Jx: (2 x 1 cell of R^2 matrix) gradient of image I2 in (row;col)
% w: (odd scalar) side length of the square around the v to evaluate C
%
% OUTPUTS
% (d x w x w) matrix of gradient evaluated at every point in (col;row)
%
% CALLEE functions
%   pointValue
%
% @author Roger Zou
% @date 6/13/15

% dimension
d = 2;

% get value of gradient at each point in cube
cJgr = pointValue(v, Jx{1}, w);
cJgc = pointValue(v, Jx{2}, w);

% compute derivative of vertex residual wrt vertex v
DfvertDv = NaN(d, w, w);
DfvertDv(1,:,:) = -cJgc;
DfvertDv(2,:,:) = -cJgr;

end
