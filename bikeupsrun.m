z=double(imread('bike2.png'))/255;
y=z(1:180,288:467);
xref=y;
y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
[k1,t1]=main_isotropic_ups(y,xref);
[k2,t2]=main_upwind_ups(y,xref);
[k3,t3]=main_cond_ups(y,xref);
[k4,t4]=main_new2_ups(y,xref);
subplot(3,2,1)
imshow(xref)
xlabel('Reference image');
subplot(3,2,2)
imshow(y)
xlabel(sprintf('Downscaled image'));
subplot(3,2,3)
imshow(k1)
xlabel(sprintf('Isotropic TV, psnr=%2.4f',t1));
subplot(3,2,4)
imshow(k2)
xlabel(sprintf('Upwind TV, psnr=%2.4f',t2));
subplot(3,2,5)
imshow(k3)
xlabel(sprintf('Condat TV, psnr=%2.4f',t3));
subplot(3,2,6)
imshow(k4)
xlabel(sprintf('New proposed model, psnr=%2.4f',t4));





