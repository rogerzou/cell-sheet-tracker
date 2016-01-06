function [ I, CC ] = drawImage( V, E, sizeI )
%DRAWIMAGE Convert from graph representation to binary image
%
% INPUTS
% V: (cell of vertices) vertex set.
% E: (cell of edges) edge set.
% sizeI: [i,j] size of image domain.
%
% OUTPUTS
% I: binary image of tracking result.
% CC: (optional) connected component from bwconncomp.
%
% @author Roger Zou
% @date 8/6/15

I = false(sizeI);

% add vertices to I
VMat = round(cell2mat(V'))';
VMat(isnan(VMat)) = [];
VInd = sub2ind(sizeI, VMat(:,2), VMat(:,1));
I(VInd) = true;

% add edges to I
for ii=1:numel(E)
    Ei = E{ii};
    if ~isempty(Ei)
        xy = Ei.curve;
        % spline interpolate for finer sampling
        n = Ei.n;
        t = 1:n;
        ts = 1:0.1:n;
        xys = spline(t,xy,ts);
        eMat = round(xys)';
        eInd = sub2ind(sizeI, eMat(:,2), eMat(:,1));
        I(eInd) = true;
    end
end

% skeletonize
I = bwmorph(I, 'thin', Inf);

% (optionally) output connected components
if nargout > 1
    CC = bwconncomp(~I, 4);
end


end
