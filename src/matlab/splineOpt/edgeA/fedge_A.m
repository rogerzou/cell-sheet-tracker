function [ f_ij, dsI, dxI ] = fedge_A( sV1, sV2, structA, structC )
%FEDGE_A Edge residual at each point for free splines.
%
% INPUTS
% sU: (struct) spline in image I1.
% sV: (struct) spline in image I2.
% structA: images struct.
%   - I1: (R^2 matrix) image I1, the previous image in a sequence.
%   - I2: (R^2 matrix) image I2, the current image in a sequence.
% structC: parameters struct.
%   - l: (odd scalar) width of spline in # of pixels.
%
% OUTPUTS
% f_ij: (double n x l) vertex residual at each point.
%
% CALLEE functions
%   splineValue
%
% @author Roger Zou
% @date 6/4/15

% get input variables
I1 = structA.I1;
J2 = structA.I2;
l = structC.l;

% compute residual
[edgeI, dsI, dxI] = splineValue(sV1, I1, l);
edgeJ = splineValue(sV2, J2, l);
f_ij = edgeI - edgeJ;

end
