% This Demo repeats the results of the follwing papers 
% Please cite them if you find them useful for your research 

% Reference:
% 1: W Lu, J Duan, Z Qiu, Z Pan, W Ryan Liu, L Bai
%    Implementation of high order variational models made easy for image processing

% 2: J Duan, Z Qiu, W Lu, G Wang, Z Pan, L Bai
%    An edge-weighted second order variational model for image decomposition

% code Writen by 
% Wenqi Lu and Jinming Duan
% contact email: j.duan@imperial.ac.uk
% March 2018

function [u,vala,valb]=TL(f,ref,s,Nbiter)
% ========================================================
% TL algorithm for denoising gray scaled images
% ========================================================
% INPUT ARGUMENTS:
% f- A noisy double graysceled image (an mXn matrix with 
% elements in [0,255])
% ref - reference grayscaled image (an mXn matrix with 
% elements in [0,255])
% s- regularization parameter
% Nbiter- number of itarations
% ========================================================
% OUTPUT ARGUMENTS:
% u- the obtained denoised image
% vala- is a 1XNbiter vector; vala(iter) is psnr(x,ref), where x is obtained denoised 
% image at step iter
% valb- is a 1XNbiter vector; valb(iter) is ssim(x,ref), where x is obtained denoised 
% image at step iter
% ========================================================
padNum=5;
f0=padarray(f,[padNum,padNum],'symmetric');
[m,n]=size(f0);
f0=double(f0);
alfa=s;theta=10;

w=zeros(m,n);
b=zeros(m,n);

u=f0;
[Y,X]=meshgrid(0:n-1,0:m-1);
G=cos(2*pi*X/m)+cos(2*pi*Y/n)-2;
vala=zeros(Nbiter,1);
valb=zeros(Nbiter,1);
for step=1:Nbiter
    %update u using FFT
    p=Bx(Fx(w-b))+By(Fy(w-b));
    g=f0+theta*p;
    u=real(ifft2(fft2(g)./(1+4*theta*G.^2)));
    
    %update w using solt thresholding 
    lap_u=Bx(Fx(u))+By(Fy(u));
    c=lap_u+b;
    q=alfa/theta;
    w=max(abs(c)-q,0).*sign(c);

    % Update Bregman parameters b
    b=c-w;
vala(step)=psnr(u(padNum+1:m-padNum,padNum+1:n-padNum)/255,ref/255);
valb(step)=ssim(u(padNum+1:m-padNum,padNum+1:n-padNum)/255,ref/255);

end
u=u(padNum+1:m-padNum,padNum+1:n-padNum)/255;
% Forward derivative operator on x with boundary condition u(:,:,1)=u(:,:,1)
function Fxu = Fx(u)
Fxu = circshift(u,[0 -1])-u;

% Forward derivative operator on y with boundary condition u(1,:,:)=u(m,:,:)
function Fyu = Fy(u)
Fyu = circshift(u,[-1 0])-u;

% Backward derivative operator on x with boundary condition Bxu(:,1)=u(:,1)
function Bxu = Bx(u)
Bxu = u - circshift(u,[0 1]);

% Backward derivative operator on y with boundary condition Bxu(1,:)=u(1,:)
function Byu = By(u)
Byu = u - circshift(u,[1 0]);









   

    