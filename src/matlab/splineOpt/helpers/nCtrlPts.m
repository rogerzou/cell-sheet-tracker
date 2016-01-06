function [ k ] = nCtrlPts( scurve, spacing )
%NCTRLPTS Computes the number of internal control points.
%
% INPUTS
% scurve: s.curve, where s is a spline.
% spacing: (scalar) spacing between each control point in pixels.
%
% OUTPUTS
% k: number of internal spline control points.
%
% CALLEE functions
%   arcLength
%
% @author Roger Zou
% @date 6/5/15

% approximate spline arc length.
len = arcLength(scurve);

% dynamically compute number of internal control points (at least one)
k = 1 + floor(len/spacing);

end
