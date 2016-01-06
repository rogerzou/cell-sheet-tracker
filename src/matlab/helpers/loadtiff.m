function [ IMatrix ] = loadtiff( fname )
%LOADTIFF Converts tiff stack file to 3D image matrix
% [ IMatrix ] = loadtiff( fname ) fname is the file path of tiff stack,
% IMatrix is the 3D image matrix.

info = imfinfo(fname);
num_images = numel(info);
IMatrix = zeros([info(1).Height, info(1).Width, num_images]);
for k = 1:num_images
    A = imread(fname, k, 'Info', info);
    IMatrix(:,:,k) = A;
end
