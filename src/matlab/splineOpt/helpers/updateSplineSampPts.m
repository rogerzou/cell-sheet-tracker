function [ s ] = updateSplineSampPts( s, interval )
%UPDATESPLINESAMPPTS recompute the sample points on spline.
%
% INPUTS
% spline: (struct) old spline.
% interval: (scalar) desired average distance between each sample point.
%
% OUTPUTS
% spline: (struct) new spline with newly spaced sample points.
%
% CALLEE functions
%   arcLength
%   splineMake
%
% @author Roger Zou
% @date 6/5/15

% approximate spline arc length
curvelength = arcLength(s.curve);

% recompute # of samples and update spline.
n = ceil(curvelength / interval);
s = splineMake(s.control, s.order, s.mult, s.open, n, false);

end
