function [x,t]=main_cond_ups(y,xref)
	
	Nbiter= 500;
	tau = 0.9/8;
	sigma = 0.9/3;
 	mu = 1;
	opD = @(x) cat(3,[diff(x,1,1);zeros(1,size(x,2))],[diff(x,1,2) zeros(size(x,1),1)]);
	opDadj = @(u) -[u(1,:,1);diff(u(:,:,1),1,1)]-[u(:,1,2) diff(u(:,:,2),1,2)];	
	prox_mu_sigma_g = @(t) t-bsxfun(@rdivide, t, max(sqrt(sum(t.^2,3))/(mu*sigma),1));
    x = prox_mu_tau_f(zeros(size(y)*2),y);
	u = zeros([size(x) 2]);
	v = zeros([size(x) 2 3]);
	tmp = opD(x);
	v(:,:,1,1) = tmp(:,:,1);
	v(:,:,2,2) = tmp(:,:,2);
	for iter = 1:Nbiter
		x = prox_mu_tau_f(x+tau*opDadj(-opD(x)+opLadj(v)-mu*u),y);
		v = prox_mu_sigma_g(v-sigma*opL(-opD(x)+opLadj(v)-mu*u));
		u = u-(-opD(x)+opLadj(v))/mu;
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
