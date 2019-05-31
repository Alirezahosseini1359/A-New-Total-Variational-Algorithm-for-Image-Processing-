function [x,vala,valb] = TGV_primal_dual(y,lambda1,lambda2,tau,Nbiter,ref)
padNum=5;
    y=padarray(y,[padNum,padNum],'symmetric');
    [m,n]=size(y);
	rho = 1.99;		% relaxation parameter, in [1,2)
	sigma = 1/tau/72; % proximal parameter
	[H,W]=size(y);
	
	opD = @(x) cat(3,[diff(x,1,1);zeros(1,W)],[diff(x,1,2) zeros(H,1)]);
	opDadj = @(u) [-u(1,:,1);-diff(u(1:end-1,:,1),1,1);u(end-1,:,1)]+...
		[-u(:,1,2) -diff(u(:,1:end-1,2),1,2) u(:,end-1,2)];	
	opJ = @(r) cat(3,...
		[r(1,:,1);diff(r(:,:,1),1,1)],[diff(r(:,:,1),1,2) zeros(H,1)],...
		[r(:,1,2) diff(r(:,:,2),1,2)],[diff(r(:,:,2),1,1);zeros(1,W)]);
	opJadj = @(u) cat(3,...
		[-diff(u(:,:,1),1,1);u(end,:,1)]-[u(:,1,2) diff(u(:,:,2),1,2)],...
		[-diff(u(:,:,3),1,2) u(:,end,3)]-[u(1,:,4);diff(u(:,:,4),1,1)]);		
	prox_tau_fx = @(x) (x+tau*y)/(1+tau);
	prox_tau_fr = @(r) r-bsxfun(@rdivide,r,max(sqrt(sum(r.^2,3))/...
		(tau*lambda1),1));
	prox_sigma_g_conj = @(u) bsxfun(@rdivide,u,max(sqrt(sum(u.^2,3))/...
		lambda2,1));
	
	x2 = y; 		% Initialization of the solution
	r2 = zeros([H,W,2]);	% Initialization of the vector field r
	u2 = zeros([H,W,4]);	% Initialization of the dual solution
		vala=zeros(Nbiter,1);
        valb=zeros(Nbiter,1);
%	fprintf('0 %f\n',lambda*sum(sum(sqrt(sum(max(opD(y),0).^2,3)))));
	for iter = 1:Nbiter
		tmp = tau*opJadj(u2);
		x = prox_tau_fx(x2-opDadj(tmp));
		r = prox_tau_fr(r2+tmp);
		u = prox_sigma_g_conj(u2+sigma*opJ(opD(2*x-x2)-(2*r-r2)));
		x2 = x2+rho*(x-x2);
		r2 = r2+rho*(r-r2);
		u2 = u2+rho*(u-u2);
        vala(iter)=psnr(x(padNum+1:m-padNum,padNum+1:n-padNum),ref);
valb(iter)=ssim(x(padNum+1:m-padNum,padNum+1:n-padNum),ref);
    end
    x=x(padNum+1:m-padNum,padNum+1:n-padNum);	
end