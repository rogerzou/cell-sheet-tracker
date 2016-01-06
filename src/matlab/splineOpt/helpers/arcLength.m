function [ len ] = arcLength( scurve )
%ARCLENGTH Computes approximate length of spline.
%
% INPUTS
% scurve: s.curve, where s is a spline.
%
% OUTPUTS
% len: approximate (arc) length of spline.
%
% CALLEE functions
%
% @author Roger Zou
% @date 6/5/15

% approximate spline arc length
curvenorm = sqrt( sum( (scurve(:,1:end-1) - scurve(:,2:end)).^2, 1) );
len = sum(curvenorm);

end
