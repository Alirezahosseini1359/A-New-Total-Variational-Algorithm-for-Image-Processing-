
% code Writen by 
% Alireza Hosseini
% contact email: hosseini.alireza@ut.ac.ir
% March 2019
function [x,t,val1,val2]=ups_tgv(y,xref,Nbiter)
% ========================================================
% TGV algorithm for upscaling gray scaled images
% ========================================================
% INPUT ARGUMENTS:
% y- A downscaled image (an mXn matrix with 
% elements in [0,1])
% ref - reference grayscaled image (an mXn matrix with 
% elements in [0,1])
% Nbiter- number of itarations
% ========================================================
% OUTPUT ARGUMENTS:
% x- the obtained denoised image
% vala- is a 1XNbiter vector; vala(iter) is psnr(x,ref), where x is obtained denoised 
% image at step iter
% valb- is a 1XNbiter vector; valb(iter) is ssim(x,ref), where x is obtained denoised 
% image at step iter
%t-psnr value at the lates iteration
% ========================================================
    lambda=1;
	alpha0=2;
    alpha1=1;
	tau = 0.9/10;
	sigma = 0.9/10;
	mu = 1;
    x = prox_mu_tau_f(zeros(size(y)*2),y);
	wbar = zeros([size(x) 2]);
    vbar = zeros([size(x) 4]);
    zbar=zeros(size(x));
    tbar=zeros([size(x) 2]);
    w=zeros([size(x) 2]);
    v=zeros([size(x) 4]);
    s=zeros(size(x));
	%fprintf('0 %f\n',sum(abs(tmp(:)))); % == sum(sum(sum(sqrt(sum(v.^2,3)))))
    val1=zeros(Nbiter,1);
        val2=zeros(Nbiter,1);

	for iter = 1:Nbiter
       [a,b,c]=opLadjco0(vbar, wbar, zbar, tbar);
       [a1,b1,c1]=opDco(x);
       x = prox_mu_tau_f(x-tau*opDadj(a1-a+mu*v,b1-b+mu*w,c1-c+mu*s),y);
       [a1,b1,c1]=opDco(x);
       a=a1-a+mu*v;
       b=b1-b+mu*w;
       c=c1-c+mu*s;
       [a1,b1,c1,d1]=opLco0(a,b,c);
    	[vbar,wbar,zbar,tbar] = prox_mu_sigma_gco(vbar+sigma*a1,wbar+sigma*b1,zbar+sigma*c1,tbar+sigma*d1 ,mu, sigma,alpha0,alpha1,lambda);
		[a,b,c] = opLadjco0(vbar,wbar,zbar,tbar);
        [d,e,f]=opDco(x);
        v = v+(d-a)/mu;
        w = w+(e-b)/mu;
        s = s+(f-c)/mu;
val1(iter)=psnr(x,xref);
val2(iter)=ssim(x,xref);

    end
  t=psnr(x,xref);
end


function xout = prox_mu_tau_f(x,y)
	z=y-(x(1:2:end,1:2:end)+x(2:2:end,1:2:end)+...
	x(1:2:end,2:2:end)+x(2:2:end,2:2:end))/4;
	xout = x;
	xout(1:2:end,1:2:end)=x(1:2:end,1:2:end)+z;
	xout(2:2:end,1:2:end)=x(2:2:end,1:2:end)+z;
	xout(1:2:end,2:2:end)=x(1:2:end,2:2:end)+z;
	xout(2:2:end,2:2:end)=x(2:2:end,2:2:end)+z;
end 
function [v,w,s]=opDco(x)
s=x;
  w=zeros([size(x) 2]);
  v=zeros([size(x) 4]);
end
function [vbar, wbar, zbar, tbar]=opLco0(v,w,s)
vbar=v;
wbar=w;
zbar=-div2(w)+s;
tbar=-div1(v)+w;
end
function x=opDadj(w,v,s)
x=s;
end
function [v,w,s]=opLadjco0(vbar, wbar, zbar, tbar)
w=wbar+Dplus(zbar)+tbar;
v=vbar-DD(tbar);
s=zbar;
end
function w=div1(v)
[m,n]=size(v(:,:,1));
	w = zeros([m n 2]);
    w(:,:,1)=Dxplus(v(:,:,1))+Dyplus(v(:,:,3));
    w(:,:,2)=Dxplus(v(:,:,4))+Dyplus(v(:,:,2));
end
function d=div2(w)
    d=Dxminus(w(:,:,1))+Dyminus(w(:,:,2));
end
function [vbar1, wbar1, sbar1, tbar1]=prox_mu_sigma_gco(vbar, wbar, sbar, tbar, mu, sigma,alpha0,alpha1,lambda)
 vbar1=vbar-bsxfun(@rdivide, vbar, max(sqrt(sum(vbar.^2,3))/(lambda*mu*sigma*alpha0),1));
 wbar1=wbar-bsxfun(@rdivide, wbar, max(sqrt(sum(wbar.^2,3))/(lambda*mu*sigma*alpha1),1));
 sbar1=sbar;
 tbar1=tbar;
end
function v=DD(w)
v(:,:,1)=-Dxminus(w(:,:,1));
v(:,:,3)=-Dyminus(w(:,:,1));
v(:,:,4)=-Dxminus(w(:,:,2));
v(:,:,2)=-Dyminus(w(:,:,2));
end
function w=Dxminus(u)
    w= [u(1:end-1,:); zeros(1,size(u,2))]-[zeros(1,size(u,2));u(1:end-1,:)];
end

function w=Dyminus(u)
    w= [u(:,1:end-1) zeros(size(u,1),1)]-[zeros(size(u,1),1) u(:,1:end-1)];
end

function w=Dplus(s)
w(:,:,1)=Dxplus(s);
w(:,:,2)=Dyplus(s);
end
function w=Dxplus(u)
    w=[u(2:end,1:end);u(end,:)]-u;
end
function w=Dyplus(u)
    w=[u(1:end,2:end) u(:,end)]-u;
end