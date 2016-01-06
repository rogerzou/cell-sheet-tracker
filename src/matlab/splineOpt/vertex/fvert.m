function [ f_i ] = fvert( u, v, structA, structC )
%FVERT Vertex residual at each point.
%
% INPUTS
% u: vertex 1.
% v: vertex 2.
% structA: images struct.
%   - I1: (R^2 matrix) image I1, the previous image in a sequence.
%   - I2: (R^2 matrix) image I2, the current image in a sequence.
% structC: parameters struct.
%   - w: (odd scalar) width of vertex window in # of pixels.
%
% OUTPUTS
% f_i: (double w x w matrix) vertex residual at each point.
%
% CALLEE functions
%   pointValue
%
% @author Roger Zou
% @date 6/13/15

% get input variables
I1 = structA.I1;
I2 = structA.I2;
w = structC.w;

% compute residual
[vertI,~] = pointValue(u, I1, w);
[vertJ,~] = pointValue(v, I2, w);
f_i = double(vertI - vertJ);

end
