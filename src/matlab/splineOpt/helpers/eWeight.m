function [ eweight ] = eWeight( structC )
%EWEIGHT Computes 2D gaussian weight matrix for edges.
% NOTE: gaussian weight supports 2.5*SD.
%
% INPUTS
% structC: parameters data structure.
%   - l: (odd scalar) spread of the spline edge in # of pixels.
%
% OUTPUTS
% vWeight: (l x 1) weight vector.
%
% CALLEE functions
%   gauss
%
% @author Roger Zou
% @date 6/2/15

% edge weights (gaussian)
halfl = floor(structC.l/2);
x = -halfl:halfl;
eweight = gauss(x, 0, halfl/2.5, true);
