function [ vweight ] = vWeight( structC )
%VWEIGHT Computes 2D gaussian weight matrix for vertices.
% NOTE: gaussian weight supports 2.5*SD.
%
% INPUTS
% structC: parameters data structure.
%   - w: (odd scalar) length of Lucas-Kanade window.
%
% OUTPUTS
% vWeight: (w x w) 2D weight matrix.
%
% CALLEE functions
%   gauss
%
% @author Roger Zou
% @date 6/2/15

% vertex weights (gaussian)
halfw = floor(structC.w/2);
x_v = -halfw:halfw;
vweight = gauss(x_v, 0, halfw/2.5, true);
vweight = vweight'*vweight;
