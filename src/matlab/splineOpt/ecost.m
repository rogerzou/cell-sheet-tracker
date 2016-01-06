function [ c ] = ecost(sV1, sV2, structA, structC, structD)
%ECOST Spline edge cost.
%
% INPUTS
% sV1: (struct) spline in image I1.
% sV2: (struct) spline in image I2.
% structA: images struct.
% structC: parameters struct.
% structD: functions struct.
%   - fedge: residual at every point along spline.
%
% OUTPUTS
% c: edge cost of spline.
%
% CALLEE functions
%   eWeight
%   integrate
%
% @author Roger Zou
% @date 7/29/15

% get input parameters
fedge = structD.fedge;

% edge weights (gaussian)
eweight = eWeight(structC);

[f_ij, dsI, dxI] = fedge(sV1, sV2, structA, structC);
costE = bsxfun(@times, (f_ij).^2, eweight);
c = integrate(costE, dsI, dxI);
