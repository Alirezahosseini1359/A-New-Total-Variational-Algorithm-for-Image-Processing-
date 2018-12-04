z=double(imread('fruits.png'))/255;
zz=rgb2gray(z);
y=zz(300:511,300:511);
xref=y;
%imwrite(xref,'bikeref.png');
y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
%figure(1);
%imshow(xref);
%colormap gray
%figure(2);
%y=y+randn(size(y))*0.02;
%imshow(y);
%colormap gray
%imwrite(y,'y.tif');
%
[k1,t1]=main_isotropic_ups(y,xref);
[k2,t2]=main_upwind_ups(y,xref);
[k3,t3]=main_cond_ups(y,xref);
[k4,t4]=main_new2_ups(y,xref);
subplot(2,3,1)
imshow(xref)
xlabel('Reference image');
subplot(2,3,2)
imshow(y)
xlabel(sprintf('Downscaled image'));
subplot(2,3,3)
imshow(k1)
xlabel(sprintf('Isotropic TV, psnr=%2.4f',t1));
subplot(2,3,4)
imshow(k2)
xlabel(sprintf('Upwind TV, psnr=%2.4f',t2));
subplot(2,3,5)
imshow(k3)
xlabel(sprintf('Condat TV, psnr=%2.4f',t3));
subplot(2,3,6)
imshow(k4)
xlabel(sprintf('New proposed model, psnr=%2.4f',t4));





