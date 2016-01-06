function [ step ] = SIFTflow( I1, I2, V )
%SIFTFLOW Compute SIFT flow to determine initial vertex positions in I2.
%
%
%


d = 2;
N = numel(V);

cellsize=3;
gridspacing=1;

sift1 = mexDenseSIFT(I1,cellsize,gridspacing);
sift2 = mexDenseSIFT(I2,cellsize,gridspacing);

SIFTflowpara.alpha=2*255;
SIFTflowpara.d=40*255;
SIFTflowpara.gamma=0.005*255;
SIFTflowpara.nlevels=4;
SIFTflowpara.wsize=2;
SIFTflowpara.topwsize=10;
SIFTflowpara.nTopIterations = 60;
SIFTflowpara.nIterations= 30;

[vx,vy,~]=SIFTflowc2f(sift1,sift2,SIFTflowpara);

Varray = cell2mat(V')';
step = zeros(N*d, 1);
step(1:2:end-1) = interp2(vx, Varray(:,1), Varray(:,2));
step(2:2:end) = interp2(vy, Varray(:,1), Varray(:,2));

step(isnan(step)) = [];

end