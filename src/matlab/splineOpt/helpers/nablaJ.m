function [ nablaj ] = nablaJ( spline, Jx, l )
%NABLAJ Gradient of J evaluated at spline points.
%
% INPUTS
% spline: (struct) spline
% Jx: 
% l: 
%
% OUTPUTS
% nablaj: 
%
% CALLEE functions
%   splineValue
%
% @author Roger Zou
% @date 6/2/15

cJgr = splineValue(spline, Jx{1}, l);
cJgc = splineValue(spline, Jx{2}, l);
nablaj = cat(3, cJgc, cJgr);
