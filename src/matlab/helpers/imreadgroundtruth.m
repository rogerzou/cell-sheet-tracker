function [ G ] = imreadgroundtruth( fpath, morph )
%IMREADGROUNDTRUTH reads the ground truth and optionally performs
%morphological operations to skeletonize.
%
% INPUTS
% fpath: (string) path of the ground truth image. Considers any pixel
% greater than 0 to be true, and 0 to be false.
% morph: (boolean) performs morphological operations to skeletonize.
%
% OUTPUTS
% G: (logical matrix) ground truth image.
%
% @author Roger Zou
% @date 2/26/15

% Validate inputs
if nargin < 2 || isempty(morph)
    morph = false;
end

% read image, convert to grayscale
if ischar(fpath)
    G = double(imread(fpath));
    if size(G, 3) == 3
        G = rgb2gray(G);
    elseif size(G, 3) == 1
    else
        disp('must be single image of grayscale or rgb');
        return
    end
else
    G = fpath;
end

% warning if more than two pixel colors in ground truth
ghist = hist(G(:));
ghist(ghist==0) = [];
if length(ghist) > 2
    warning('IMREADGROUNDTRUTH: image from path: %s is not a binary image', fpath);
end

% convert to logical type
G(G>=1) = 1;
G = logical(G);

% if 'morph'==true, then thin
if morph
    G = bwmorph(G, 'thin', Inf);
end
