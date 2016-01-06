function spline = updateSplineEndPts_A(spline, a, b)
%UPDATESPLINEENDPTS_A Compute updated spline s with new endpoints.
%
% INPUT
% spline: (struct) spline.
% a: (2x1 vector) of endpoint a.
% b: (2x1 vector) of endpoint b.
%
% OUTPUT
% spline: (struct) updated spline with new control points.
%
% CALLEE functions
%   splineEvalEven
%
% @author Roger Zou
% @date 8/15/15

% get old spline endpoints
olda = spline.control(:,1);
oldb = spline.control(:,end);

% compute rotation matrix of axis joining a and b
oldeu = getE(olda,oldb);
neweu = getE(a, b);
theta = atan2(det([oldeu, neweu]), dot(oldeu, neweu));
R = [cos(theta),-sin(theta); sin(theta), cos(theta)];

% rotate interior control points (centered around a) accordingly
oldintctrl = bsxfun(@minus, spline.control(:,2:end-1), olda);
newintctrl = bsxfun(@plus, R*oldintctrl, a);
spline.control = [a,newintctrl,b];
spline = splineEvalEven(spline, true, true, false);

end
