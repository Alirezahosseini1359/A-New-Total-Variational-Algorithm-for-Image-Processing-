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

function [u,vala,valb]=infcon(f,ref,s,Nbiter)
% ========================================================
% infcon algorithm for denoising gray scaled images
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
f0=double(f0);
[m,n]=size(f0);

%%%%%%%
u2=f0;
%%%%%%%

b1=zeros(m,n);
b2=zeros(m,n);
w1=zeros(m,n);
w2=zeros(m,n);
s1=zeros(m,n);
s2=zeros(m,n);
s3=zeros(m,n);
d1=zeros(m,n);
d2=zeros(m,n);
d3=zeros(m,n);

alfa=s;    %(lamda)
beta=2*s;    %(mu)
theta1=5;
theta2=5;

[Y,X]=meshgrid(0:n-1,0:m-1);
G=cos(2*pi*X/m)+cos(2*pi*Y/n)-2;
vala=zeros(Nbiter,1);
valb=zeros(Nbiter,1);
for step=1:Nbiter  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update u using FFT
    div_w_b=Bx(w1-b1)+By(w2-b2);
    g=f0-u2-theta1*div_w_b;
    u1=real(ifft2(fft2(g)./(1-2*theta1*G)));
    
    % update w using solt thresholding
    c1=Fx(u1)+b1;
    c2=Fy(u1)+b2;
    abs_c=sqrt(c1.^2+c2.^2+eps);
    w1=max(abs_c-alfa/theta1,0).*c1./abs_c;
    w2=max(abs_c-alfa/theta1,0).*c2./abs_c;
    
    % update Bregman iterative parameters
    b1=c1-w1;
    b2=c2-w2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    div_w_b=Bx(Fx(s1-d1))+By(Fy(s2-d2))+2*By(Bx(s3-d3));
    g=f0-u1+theta2*div_w_b;
    u2=real(ifftn(fftn(g)./(1+4*theta2*G.^2)));
    
    c1=Bx(Fx(u2))+d1;
    c2=By(Fy(u2))+d2;
    c3=Fy(Fx(u2))+d3;
    abs_c=sqrt(c1.^2+c2.^2+2*c3.^2+eps);
    s1=max(abs_c-beta/theta2,0).*c1./abs_c;
    s2=max(abs_c-beta/theta2,0).*c2./abs_c;
    s3=max(abs_c-beta/theta2,0).*c3./abs_c;
    
    % update Bregman iterative parameters
    d1=c1-s1;
    d2=c2-s2;
    d3=c3-s3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  v=u1+u2;
v=v(padNum+1:m-padNum,padNum+1:n-padNum);  
vala(step)=psnr(v/255,ref/255);
valb(step)=ssim(v/255,ref/255);
end
u=u1+u2;
u=u(padNum+1:m-padNum,padNum+1:n-padNum)/255;
% Forward derivative operator on x with periodic boundary condition
function Fxu = Fx(u)
Fxu = circshift(u,[0 -1])-u;

% Forward derivative operator on y with periodic boundary condition
function Fyu = Fy(u)
Fyu = circshift(u,[-1 0])-u;

% Backward derivative operator on x with periodic boundary condition
function Bxu = Bx(u)
Bxu = u-circshift(u,[0 1]);

% Backward derivative operator on y with periodic boundary condition
function Byu = By(u)
Byu = u-circshift(u,[1 0]);
