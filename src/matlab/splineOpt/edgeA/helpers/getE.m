function [ eu, ev, DeuDabT, DevDabT ] = getE( a, b )
%GETE get spline coordinate vectors {eu,ev} and derivatives wrt endpoints
%
% INPUT
% a: (d x 1) endpoint 1.
% b: (d x 1) endpoint 2.
%
% OUTPUT
% eu: (d x 1) unit vector along spline.
% ev: (d x 1) unit vector rotated 90 deg counter-clockwise from eu.
% DeuDvT: (d x 2d) jacobian matrix of derivative of eu wrt both endpoints.
% DevDvT: (d x 2d) jacobian matrix of derivative of ev wrt both endpoints.
%
% @author Roger Zou
% @date 6/9/15

d = length(a);

% get difference and norm
diffS = b - a;
normS = norm(b - a);

% get eu and ev
R = [0,-1; 1,0];
eu = diffS / normS;
ev = R * eu;

% get jacobians
if nargout > 2
    DeuDaT = diffS * diffS' / normS^3 - eye(d) / normS;
    DeuDabT = { DeuDaT, - DeuDaT };
    DevDabT = { R * DeuDaT, - R * DeuDaT };
end
