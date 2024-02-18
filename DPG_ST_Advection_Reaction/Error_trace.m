function I = Error_trace(uexact,uhat,x0,x1,t1,order_s)
%Change this routine using quadratures

h=x1-x0;
xm=(x1+x0)/2;
LT=@(x) 2*(x-xm)/h;
Q0=@(x) 1;
Q1=@(x) x;
Q2=@(x) (3*x.^2-1)/2;
Q3=@(x) (5*x.^3-3*x)/2;

if order_s==0
    u_num=@(x) (uexact(x,t1)-uhat.*Q0(x)).^2;
    I=integral(u_num,x0,x1);
elseif order_s==1
    u_num=@(x) (uexact(x,t1)-(uhat(1).*Q0(x)+uhat(2).*Q1(LT(x)))).^2;
    I=integral(u_num,x0,x1);
elseif order_s==2
    u_num=@(x) (uexact(x,t1)-(uhat(1).*Q0(x)+uhat(2).*Q1(LT(x))+uhat(3).*Q2(LT(x)))).^2;
    I=integral(u_num,x0,x1);
elseif order_s==3
    u_num=@(x) (uexact(x,t1)-(uhat(1).*Q0(x)+uhat(2).*Q1(LT(x))+uhat(3).*Q2(LT(x))+uhat(4).*Q3(LT(x)))).^2;
    I=integral(u_num,x0,x1);
end

end

