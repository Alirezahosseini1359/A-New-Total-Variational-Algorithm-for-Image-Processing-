%In this code, denoising test image cameraman is solved
%with 3 different variational methods:
%1-proposed 2-TGV(split Bregman) 3-TGV (primal dual) 
%================================================================================
close all; clc 
ref=double(imread('cameraman.tif'))/255;% original image 
ref=ref(100:300,100:300); % a part of image is used in comparison
f=ref+randn(size(ref))*0.18; %noisy image with additive white Gaussian noise of standard deviation 0.18.
Nbiter=500;
[x11,vala1,valb1] = TGVdenoising(f,.15,.3,0.99/10,Nbiter,ref);
[u1,vala3,valb3]=TGV(f,ref,0.15,Nbiter);
[x,vala2,valb2]=proposed(f,ref,0.08,Nbiter);
xx=1:500;
figure
plot(xx,vala1, xx,vala3, xx,vala2,'LineWidth', 1)
legend('TGV_primal_dual', 'TGV_split_Bregman', 'Proposed','Location', 'SouthEast')
title('Comparison of PSNR')
figure
plot(xx,valb1, xx,valb3, xx,valb2,'LineWidth', 1)
legend('TGV-primal-dual', 'TGV-split-Bregman', 'Proposed-primal-dual','Location', 'SouthEast')
title('Comparison of SSIM')
