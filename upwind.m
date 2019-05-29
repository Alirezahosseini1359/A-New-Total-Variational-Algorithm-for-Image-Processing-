function [x,vala,valb]=upwind(y,ref,lambda,Nbiter) 
% ========================================================
% upwind algorithm for denoising gray scaled images
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
	tau = 0.01;
	sigma = 1/tau/16;
	rho = 1.9;		
	opD = @(x) cat(3,[-diff(x,1,1);zeros(1,size(x,2))],[-diff(x,1,2) zeros(size(x,1),1)],...
		[zeros(1,size(x,2));diff(x,1,1)],[zeros(size(x,1),1) diff(x,1,2)]);
	opDadj = @(u) [u(1,:,1);diff(u(:,:,1),1,1)]+[u(:,1,2) diff(u(:,:,2),1,2)]+...
		[-diff(u(:,:,3),1,1);u(end,:,3)]+[-diff(u(:,:,4),1,2) u(:,end,4)];	
	prox_tau_f = @(x) (x+tau*y)/(1+tau);
	
	x = y;
	u = zeros([size(x) 4]);
    vala=zeros(Nbiter,1);
        valb=zeros(Nbiter,1);
	for iter = 1:Nbiter
		xnew = prox_tau_f(x-tau*opDadj(u));
		unew = prox_sigma_g_conj(u+sigma*opD(2*xnew-x),lambda);
		x = xnew+(rho-1)*(xnew-x);
		u = unew+(rho-1)*(unew-u);
vala(iter)=psnr(x(padNum+1:m-padNum,padNum+1:n-padNum),ref);
valb(iter)=ssim(x(padNum+1:m-padNum,padNum+1:n-padNum),ref);
	end
x=x(padNum+1:m-padNum,padNum+1:n-padNum);	
end

function unew = prox_sigma_g_conj(u,lambda) 
	u = max(u,0);
	unew = bsxfun(@rdivide, u, max(sqrt(sum(u.^2,3))/lambda,1));
end