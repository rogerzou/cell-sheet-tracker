function [ s ] = estep_A( s, step )
%ESTEP_A Step in the internal spline control points.
%
% INPUTS
% s: (struct) old spline.
% step: (r x 1) vector of step values for internal independent variables.
%
% OUTPUTS
% s: (struct) new spline after taking step.
%
% CALLEE functions
%   splineEvalEven
%
% @author Roger Zou
% @date 6/3/15

% take newton step
muq = s.control(:, 2:end-1);
[~, k] = size(muq);
muq = reshape(muq(:)+step, s.d, k);
control = [s.control(:,1), muq, s.control(:,end)];

% update spline
s.control = control;
s = splineEvalEven(s, true, true, false);

end
