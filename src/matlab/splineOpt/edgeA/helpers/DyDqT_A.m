function [ DyDqT ] = DyDqT_A( spline, l )
%DYDQT_A  Derivative of points along the spline wrt control points
%\frac{\partial \by}{\bq^T}
%
% INPUT
% spline: (struct spline
% l: (odd scalar) spread of the spline edge in # of pixels
%
% OUTPUT
% DyDqT: (d x m x n x l) matrix of the derivative DyDq at each point along
% spline.
%
% @author Roger Zou
% @date 6/3/15

% dx
halfl = floor(l/2);
x = -halfl:halfl;
% get sizes - d  = dimension, m = 2(k+2)
[d, n] = size(spline.normal);
[~, m] = size(spline.control);
m = m * d;

DyDqT = NaN(d, m, n, l);
for i=1:n
    for j=1:l
        indexrange = (2*i)-1 : (2*i)-1+(d-1);
        DyDqT(:,:,i,j) = spline.P(indexrange,:) + spline.N(indexrange,:) * x(j);
    end
end


end
