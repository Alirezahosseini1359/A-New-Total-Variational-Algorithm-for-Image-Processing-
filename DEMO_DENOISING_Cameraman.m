%In this code, denoising test image cameraman is solved
%with 8 different variational methods:
%1- TC 2-proposed 3-TGV 4-TL 5-Condat  6-Isotropic TV 7-Upwind 8-Infcon
%================================================================================
close all; clc 
ref=double(imread('cameraman.tif'))/255;% original image 
ref=ref(100:300,100:300); % a part of image is used in comparison
sTC=10;%regularization parameter for TC
sNew=0.08;%regularization parameter for proposed model
sTGV=0.16; %regularization parameter for TGV
sTL=50; %regularization parameter for TL
sCondat=0.14; %regularization parameter for Condat
sTV=0.15; %regularization parameter for Isotropic TV
sUpwind=0.19; %regularization parameter for Upwind TV
sInf=40; %regularization parameter for Infcon
f=ref+randn(size(ref))*0.18; %noisy image with additive white Gaussian noise of standard deviation 0.18.
ff=f*255; 
Nbiter=500;%number of iterations
[u,vala1,valb1]=TC(ff,ref*255,sTC,Nbiter);
[x,vala2,valb2]=proposed(f,ref,sNew,Nbiter);
[u1,vala3,valb3]=TGV(f,ref,sTGV,Nbiter);
[u2,vala4,valb4]=TL(ff,ref*255,sTL,Nbiter);
[x1,vala5,valb5]=condat(f,ref,sCondat,Nbiter);
[u3,vala6,valb6]=TV(f,ref,sTV,Nbiter);
[x3,vala7,valb7]=upwind(f,ref,sUpwind,Nbiter);
[u4,vala8,valb8]=infcon(ff,ref*255,sInf,Nbiter);
xx=1:Nbiter;
figure
plot(xx,vala1, xx,vala2, xx,vala3, xx,vala4, xx,vala5, xx,vala6, xx,vala7,xx,vala8,'LineWidth', 1)
legend('TC', 'Proposed', 'TGV', 'TL', 'Condat TV', 'Isotropic', 'Upwind', 'INFCON','Location', 'SouthEast')
title('Comparison of PSNR')

figure

plot(xx,valb1, xx,valb2, xx,valb3, xx,valb4, xx,valb5,xx,valb6,xx,valb7, xx,valb8,'LineWidth', 1)
legend('TC', 'Proposed', 'TGV', 'TL', 'Condat TV', 'Isotropic TV', 'Upwind TV', 'INFCON', 'Location', 'SouthEast')
title('Comparison of SSIM') 
figure
subplot(3,4,1)
imshow(ref)
xlabel('Reference image');
subplot(3,4,2)
imshow(f)
xlabel(sprintf('Noisy image'));
subplot(3,4,3)
imshow(u)
xlabel(sprintf('TC'));
subplot(3,4,4)
imshow(x)
xlabel(sprintf('proposed'));
subplot(3,4,5)
imshow(u1)
xlabel(sprintf('TGV'));
subplot(3,4,6)
imshow(u2)
xlabel(sprintf('TL'));
subplot(3,4,7)
imshow(x1)
xlabel(sprintf('Condat'));
subplot(3,4,8)
imshow(u3)
xlabel(sprintf('TV'));
subplot(3,4,9)
imshow(x3)
xlabel(sprintf('Upwind'));
subplot(3,4,10)
imshow(u4)
xlabel(sprintf('Infcon'));
T1=[vala1(end) vala2(end) vala3(end) vala4(end) vala5(end) vala6(end) vala7(end) vala8(end)];
T2=[valb1(end) valb2(end) valb3(end) valb4(end) valb5(end) valb6(end) valb7(end) valb8(end)];
disp(' PSNR values') 
disp(' TC     Proposed     TGV      TL    Condat     TV     Upwind     Infcon') 
disp('----------------------------------------------------------------') 
fprintf(' %2.4f  %2.4f  %2.4f  %2.4f  %2.4f  %2.4f  %2.4f' ,T1(1),T1(2),T1(3),T1(4),T1(5),T1(6),T1(7),T1(8)); 
fprintf('\n');
fprintf('\n');
disp(' SSIM values') 
disp(' TC    Proposed    TGV     TL   Condat    TV    Upwind    Infcon') 
disp('----------------------------------------------------------------') 
fprintf(' %0.4f  %0.4f  %0.4f  %0.4f  %0.4f  %0.4f  %0.4f' ,T2(1),T2(2),T2(3),T2(4),T2(5),T2(6),T2(7),T2(8));

