function [xnew,t]=main_isotropic_ups(y,xref)
	Nbiter= 500;
	tau = 1/8;
	sigma = 0.1;
	rho = 1;			
	opD = @(x) cat(3,[diff(x,1,1);zeros(1,size(x,2))],[diff(x,1,2) zeros(size(x,1),1)]);
	opDadj = @(u) [zeros(1,size(u,2)); u(1:end-1,:,1)]-[u(1:end-1,:,1);zeros(1,size(u,2))]+...
        [zeros(size(u,1),1) u(:,1:end-1,2)]-[u(:,1:end-1,2) zeros(size(u,1),1)];
	prox_sigma_g_conj = @(u) bsxfun(@rdivide, u, max(sqrt(sum(u.^2,3)),1));
	x = prox_tau_f(zeros(size(y)*2),y);
	u = zeros([size(x) 2]);
	for iter = 1:Nbiter
		xnew = prox_tau_f(x-tau*opDadj(u),y);
		unew = prox_sigma_g_conj(u+sigma*opD(2*xnew-x));
		x = xnew+(rho-1)*(xnew-x);
		u = unew+(rho-1)*(unew-u);
    end
t=psnr(xnew,xref);
end


function xout = prox_tau_f(x,y)
	z=y-(x(1:2:end,1:2:end)+x(2:2:end,1:2:end)+...
	x(1:2:end,2:2:end)+x(2:2:end,2:2:end))/4;
	xout = x;
	xout(1:2:end,1:2:end)=x(1:2:end,1:2:end)+z;
	xout(2:2:end,1:2:end)=x(2:2:end,1:2:end)+z;
	xout(1:2:end,2:2:end)=x(1:2:end,2:2:end)+z;
	xout(2:2:end,2:2:end)=x(2:2:end,2:2:end)+z;
end 

