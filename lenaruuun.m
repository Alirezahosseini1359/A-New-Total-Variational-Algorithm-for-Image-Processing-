zz=double(imread('lena_gray_512.png'))/255;
z=zz(200:300,200:400);
yy=z+randn(size(z))*0.1;
[t1,t2,t3,t4,k1,k2,k3,k4,a1,a2,a3,a4]=main_lambda(yy,z);
subplot(2,3,1)
imshow(z)
xlabel('Reference image');
subplot(2,3,2)
imshow(yy)
xlabel(sprintf('Noisy image'));
subplot(2,3,3)
imshow(k1)
xlabel(sprintf('Isotropic TV, lambda=%1.2f, psnr=%2.4f',a1,t1));
subplot(2,3,4)
imshow(k2)
xlabel(sprintf('Upwind TV, lambda=%1.2f, psnr=%2.4f',a2, t2));
subplot(2,3,5)
imshow(k3)
xlabel(sprintf('Condat TV, lambda=%1.2f, psnr=%2.4f',a3, t3));
subplot(2,3,6)
imshow(k4)
xlabel(sprintf('New proposed model, lambda=%1.2f, psnr=%2.4f',a4,t4));





