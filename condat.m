function [x,vala,valb]=condat(y,ref,lambda,Nbiter) 
% ========================================================
% Condat's algorithm for denoising gray scaled images
% ========================================================
% INPUT ARGUMENTS:
% y- A noisy double graysceled image (an mXn matrix with 
% elements in [0,1])
% ref - reference grayscaled image (an mXn matrix with 
% elements in [0,1])
% lambda- regularization parameter
% Nbiter- number of itarations
% ========================================================
% OUTPUT ARGUMENTS:
% x- the obtained denoised image
% vala- is a 1XNbiter vector; vala(iter) is psnr(x,ref), where x is obtained denoised 
% image at step iter
% valb- is a 1XNbiter vector; valb(iter) is ssim(x,ref), where x is obtained denoised 
% image at step iter
% ========================================================

	padNum=5;
    y=padarray(y,[padNum,padNum],'symmetric');
    [m,n]=size(y);
	tau = 0.99/8;
	sigma = 0.99/3;
	mu = 1; 
	% 			
	opD = @(x) cat(3,[diff(x,1,1);zeros(1,size(x,2))],[diff(x,1,2) zeros(size(x,1),1)]);
	opDadj = @(u) -[u(1,:,1);diff(u(:,:,1),1,1)]-[u(:,1,2) diff(u(:,:,2),1,2)];	
	prox_mu_tau_f = @(x) (x+mu*tau*y)/(1+mu*tau);
	prox_mu_sigma_g = @(t) t-bsxfun(@rdivide, t, max(sqrt(sum(t.^2,3))/(lambda*mu*sigma),1));
	
	x = y;
	u = zeros([size(x) 2]);
	v = zeros([size(y) 2 3]);
	tmp = opD(x);
	v(:,:,1,1) = tmp(:,:,1);
	v(:,:,2,2) = tmp(:,:,2);
    vala=zeros(Nbiter,1);
        valb=zeros(Nbiter,1);
	for iter = 1:Nbiter
		x = prox_mu_tau_f(x+tau*opDadj(-opD(x)+opLadj(v)-mu*u));
		v = prox_mu_sigma_g(v-sigma*opL(-opD(x)+opLadj(v)-mu*u));
		u = u-(-opD(x)+opLadj(v))/mu;
vala(iter)=psnr(x(padNum+1:m-padNum,padNum+1:n-padNum),ref);
valb(iter)=ssim(x(padNum+1:m-padNum,padNum+1:n-padNum),ref);
    end
x=x(padNum+1:m-padNum,padNum+1:n-padNum);
end

function t = opL(u)
	[height,width,d]=size(u);
	t=zeros(height,width,2,3);
	t(:,:,1,1)=u(:,:,1); 
	t(1:end-1,2:end,2,1)=(u(2:end,1:end-1,2)+u(1:end-1,1:end-1,2)+u(2:end,2:end,2)+u(1:end-1,2:end,2))/4; 			t(1:end-1,1,2,1)=(u(1:end-1,1,2)+u(2:end,1,2))/4; 
	t(:,:,2,2)=u(:,:,2);
	t(2:end,1:end-1,1,2)=(u(2:end,1:end-1,1)+u(1:end-1,1:end-1,1)+u(2:end,2:end,1)+u(1:end-1,2:end,1))/4; 
	t(1,1:end-1,1,2)=(u(1,1:end-1,1)+u(1,2:end,1))/4; 	
	t(2:end,:,1,3) = (u(2:end,:,1)+u(1:end-1,:,1))/2; 
	t(1,:,1,3) = u(1,:,1)/2;
	t(:,2:end,2,3) = (u(:,2:end,2)+u(:,1:end-1,2))/2; 
	t(:,1,2,3) = u(:,1,2)/2;
end

function u = opLadj(t)
	[height,width,d,c]=size(t);
	u=zeros(height,width,2);
	u(1:end-1,2:end,1)=t(1:end-1,2:end,1,1)+(t(1:end-1,2:end,1,2)+t(1:end-1,1:end-1,1,2)+...
	t(2:end,2:end,1,2)+t(2:end,1:end-1,1,2))/4+(t(1:end-1,2:end,1,3)+t(2:end,2:end,1,3))/2;
	u(1:end-1,1,1)=t(1:end-1,1,1,1)+(t(1:end-1,1,1,2)+t(2:end,1,1,2))/4+...
	(t(1:end-1,1,1,3)+t(2:end,1,1,3))/2;
	u(2:end,1:end-1,2)=t(2:end,1:end-1,2,2)+(t(2:end,1:end-1,2,1)+t(1:end-1,1:end-1,2,1)+...
	t(2:end,2:end,2,1)+t(1:end-1,2:end,2,1))/4+(t(2:end,1:end-1,2,3)+t(2:end,2:end,2,3))/2;
	u(1,1:end-1,2)=t(1,1:end-1,2,2)+(t(1,1:end-1,2,1)+t(1,2:end,2,1))/4+...
	(t(1,1:end-1,2,3)+t(1,2:end,2,3))/2;
end
