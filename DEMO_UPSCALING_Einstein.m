%In this code, upscaling test image Einstein is solved
%with 5 different variational methods:
%1- isotropic tv 2-upwind 3-condat 4-proposed 5-tgv 
%================================================================================
close all; clc 
Nbiter=300;
zz=double(imread('ensh.tif'))/255;
y=zz(1:672,:);
xref=y;
xref=(xref(1:2:end,1:2:end)+xref(2:2:end,1:2:end)+xref(1:2:end,2:2:end)+xref(2:2:end,2:2:end))/4;
%imwrite(xref,'bikeref.png');
y=(xref(1:2:end,1:2:end)+xref(2:2:end,1:2:end)+xref(1:2:end,2:2:end)+xref(2:2:end,2:2:end))/4;
[k1,t1,val1,wal1]=ups_isotropic(y,xref,Nbiter);
[k2,t2,val2,wal2]=ups_upwind(y,xref,Nbiter);
[k3,t3,val3,wal3]=ups_condat(y,xref,Nbiter);
[k4,t4,val4,wal4]=ups_proposed(y,xref,Nbiter);
[k5,t5,val5,wal5]=ups_tgv(y,xref,Nbiter);
xx=1:Nbiter;
figure
plot(xx,val1, xx,val2, xx,val3, xx,val4, xx,val5,'LineWidth', 1)
legend('TV', 'Upwind TV', 'Condat TV', 'New TV', 'TGV', 'Location', 'SouthEast')
title('Comparison of PSNR')
figure
plot(xx,wal1, xx,wal2, xx,wal3, xx,wal4, xx,wal5,'LineWidth', 1)
legend('TV', 'Upwind TV', 'Condat TV', 'New TV', 'TGV', 'Location', 'SouthEast')
title('Comparison of SSIM')

figure
subplot(3,3,1)
imshow(xref)
xlabel('Reference image');
subplot(3,3,2)
imshow(y)
xlabel(sprintf('Downscaled image'));
subplot(3,3,3)
imshow(k1)
xlabel(sprintf('Isotropic'));
subplot(3,3,4)
imshow(k2)
xlabel(sprintf('Upwind'));
subplot(3,3,5)
imshow(k3)
xlabel(sprintf('Condat'));
subplot(3,3,6)
imshow(k4)
xlabel(sprintf('Proposed'));
subplot(3,3,7)
imshow(k5)
xlabel(sprintf('TGV'));
T1=[val1(end) val2(end) val3(end) val4(end) val5(end)];
T2=[wal1(end) wal2(end) wal3(end) wal4(end) wal5(end)];
disp(' PSNR values') 
disp(' Isotropic   upwind   Condat  Proposed  TGV') 
disp('----------------------------------------------------------------') 
fprintf(' %2.4f  %2.4f  %2.4f  %2.4f  %2.4f ' ,T1(1),T1(2),T1(3),T1(4),T1(5)); 
fprintf('\n');
fprintf('\n');
disp(' SSIM values') 
disp(' Isotropic  upwind  Condat   Proposed TGV') 
disp('----------------------------------------------------------------') 
fprintf(' %0.4f  %0.4f  %0.4f  %0.4f  %0.4f' ,T2(1),T2(2),T2(3),T2(4),T2(5));

