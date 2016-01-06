function [ Values, dr, dc, A] = splineValue( spline, C, l )
%SPLINEVALUE Evaluates the value of the spline at cost function
% Evaluates the cost function at spline coordinates, with spread determined
% by a hypothetical gaussian. Note that the output is NOT WEIGHTED in any
% manner.
%
% INPUTS
% spline: (struct) a spline on cost function C.
% C: (R^2 matrix) image cost matrix function.
% l: (odd scalar) spread of the spline edge in # of pixels
%
% OUTPUTS
% Values: (n x l) matrix of the cost at each pixel considered for the
% overall spline cost, UNWEIGHTED.
% dr: (vector) differential with respect to row of Values (ds of spline)
% dc: (vector) differential with respect to column of Values (dx of spline)
% A: (n x l x d) coordinate matrix of the evaluated points in C
%
% @author Roger Zou
% @date 8/18/15

% convert types
C = double(C);

% get sizes
halfl = floor(l/2);
[r,c] = size(C);

% get spline normals and curve samples
nmal = spline.normal';
curve = spline.curve';

% determine the coordinates of n points lying on each of the m needles
llen = -halfl:halfl;
xlen = nmal(:,1) * llen;
ylen = nmal(:,2) * llen;
Ap1 = bsxfun(@plus, curve(:,1), xlen);
Ap2 = bsxfun(@plus, curve(:,2), ylen);
A = cat(3, Ap1, Ap2);

% check for out of bounds coordinates
outofbounds = (A(:,:,1) < 1 | A(:,:,2) < 1 | A(:,:,1) > c | A(:,:,2) > r);
if sum(outofbounds(:)) > 0
    warning('SPLINEVALUE: out of bounds! These values will be set to Inf.');
end

% compute costs at the coordinates of all spline needles (bilinear interp)
Vq = interp2(C, Ap1(:), Ap2(:));
Vq(isnan(Vq)) = Inf;
Values = reshape(Vq, size(Ap1));

% compute dr (differential wrt row of Values);
% dr is the vector of distances between consecutive samples on the spline 
dx = diff(Ap1(:, halfl+1));
dy = diff(Ap2(:, halfl+1));
dr = sqrt(dx .^ 2 + dy .^ 2);

% compute dc (differential wrt column of Values)
dx = diff(Ap1(1, :));
dy = diff(Ap2(2, :));
dc = sqrt(dx .^ 2 + dy .^ 2);

end
