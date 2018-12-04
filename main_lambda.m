function [t1,t2,t3,t4,x01,x02,x03,x04,a1,a2,a3,a4]=main_lambda(yy,z)
s=0.05;
x01=main_isotropic(yy,s);
d=x01;
s=s+0.01;
x=main_isotropic(yy,s);
dd=x;
if psnr(x,z)> psnr(x01,z)
       while psnr(x,z)> psnr(x01,z)
           x01=x;
           s=s+.01;
           x=main_isotropic(yy,s);
       end
       a1=s-0.01;
else

x01=d;
x=dd;
if psnr(x,z)< psnr(x01,z)
    s=s-0.01;
x=main_isotropic(yy,s);
       while psnr(x,z)> psnr(x01,z)
           x01=x;
           s=s-.02;
           x=main_isotropic(yy,s);
       end
       a1=s+0.01;
end
end
t1=psnr(x01,z);
s=0.05;
x02=main_upwind(yy,s);
s=s+0.01;
x=main_upwind(yy,s);
d=x02;
dd=x;
if psnr(x,z)> psnr(x02,z)
       while psnr(x,z)> psnr(x02,z)
           x02=x;
           s=s+.01;
           x=main_upwind(yy,s);
       end
       a2=s-0.01;
else

x02=d;
x=dd;
if psnr(x,z)< psnr(x02,z)
    s=s-0.01;
x=main_upwind(yy,s);
       while psnr(x,z)> psnr(x02,z)
           x02=x;
           s=s-.02;
           x=main_upwind(yy,s);
       end
       a2=s+0.01;
end
end
t2=psnr(x02,z);

s=0.05;
x03=main_condat(yy,s);
s=s+0.01;
x=main_condat(yy,s);
d=x03;
dd=x;
if psnr(x,z)> psnr(x03,z)
       while psnr(x,z)> psnr(x03,z)
           x03=x;
           s=s+.01;
           x=main_condat(yy,s);
       end
       a3=s-0.01;
else

x03=d;
x=dd;
if psnr(x,z)< psnr(x03,z)
    s=s-0.01;
x=main_condat(yy,s);
       while psnr(x,z)> psnr(x03,z)
           x03=x;
           s=s-.02;
           x=main_condat(yy,s);
       end
       a3=s+0.01;
end
end
t3=psnr(x03,z);
s=0.05;
x04=main_newalgorithm2(yy,s);
s=s+0.01;
x=main_newalgorithm2(yy,s);
d=x04;
dd=x;
if psnr(x,z)> psnr(x04,z)
       while psnr(x,z)> psnr(x04,z)
           x04=x;
           s=s+.01;
           x=main_newalgorithm2(yy,s);
       end
       a4=s-0.01;
else

x04=d;
x=dd;
if psnr(x,z)< psnr(x04,z)
    s=s-0.02;
x=main_newalgorithm2(yy,s);
       while psnr(x,z)> psnr(x04,z)
           x04=x;
           s=s-.01;
           x=main_newalgorithm2(yy,s);
       end
       a4=s+0.01;
end
end
t4=psnr(x04,z);
