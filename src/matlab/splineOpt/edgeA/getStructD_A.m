function [ structD ] = getStructD_A
%GETSTRUCTD_A Get structD for edgeA.
%
% @author Roger Zou
% @date 8/13/15

structD = struct('fvert', @fvert, 'DfvertDv', @DfvertDv, ...
        'fedge', @fedge_A, 'DfedgeDs', @DfedgeDs_A, ...
        'DfedgeDv', @DfedgeDv_A, 'estep', @estep_A, ...
        'optInitSpline', @optInitSpline_A, 'updateSplineEndPts', @updateSplineEndPts_A );
