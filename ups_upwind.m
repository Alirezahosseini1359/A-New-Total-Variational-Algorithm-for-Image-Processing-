function [xnew,t,val1,val2]=ups_upwind(y,xref,Nbiter)
		tau = 0.02;
	sigma = 1/tau/16;
	rho = 1.9;	
	opD = @(x) cat(3,[-diff(x,1,1);zeros(1,size(x,2))],[-diff(x,1,2) zeros(size(x,1),1)],...
		[zeros(1,size(x,2));diff(x,1,1)],[zeros(size(x,1),1) diff(x,1,2)]);
	opDadj = @(u) [u(1,:,1);diff(u(:,:,1),1,1)]+[u(:,1,2) diff(u(:,:,2),1,2)]+...
		[-diff(u(:,:,3),1,1);u(end,:,3)]+[-diff(u(:,:,4),1,2) u(:,end,4)];	
	
	x = prox_tau_f(zeros(size(y)*2),y);
	u = zeros([size(x) 4]);
    val1=zeros(1,Nbiter);
    val2=zeros(1,Nbiter);
	for iter = 1:Nbiter
		xnew = prox_tau_f(x-tau*opDadj(u),y);
		unew = prox_sigma_g_conj(u+sigma*opD(2*xnew-x));
		x = xnew+(rho-1)*(xnew-x);
 		u = unew+(rho-1)*(unew-u);
val1(iter)=psnr(xnew,xref);
val2(iter)=ssim(xnew,xref);
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

function unew = prox_sigma_g_conj(u) 
	u = max(u,0);
	unew = bsxfun(@rdivide, u, max(sqrt(sum(u.^2,3)),1));
end