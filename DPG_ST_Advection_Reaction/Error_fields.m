function I = Error_fields(uexact,u0,u1,u2,u3,x0,x1,t0,t1,order_t,order_s)
%Change this routine using quadratures
P0=@(t) 1;
P1=@(t) (t-t0)/(t1-t0);
P2=@(t) ((t-t0)/(t1-t0)).^2;
P3=@(t) ((t-t0)/(t1-t0)).^3;

h=x1-x0;
xm=(x1+x0)/2;
LT=@(x) 2*(x-xm)/h;
Q0=@(x) 1;
Q1=@(x) x;
Q2=@(x) (3*x.^2-1)/2;
Q3=@(x) (5*x.^3-3*x)/2;

if order_t==1
    if order_s==0
        u_num=@(x,t) (uexact(x,t)-u0.*P0(t).*Q0(x)).^2;
        %I=integral2(u_num,x0,x1,t0,t1);
        I=quad2d(u_num,x0,x1,t0,t1);
    elseif order_s==1
        u_num=@(x,t) (uexact(x,t)-(u0(1).*P0(t).*Q0(x)+u0(2).*P0(t).*Q1(LT(x)))).^2;
        %I=integral2(u_num,x0,x1,t0,t1);
        I=quad2d(u_num,x0,x1,t0,t1);
    elseif order_s==2
        u_num=@(x,t) (uexact(x,t)-(u0(1).*P0(t).*Q0(x)+u0(2).*P0(t).*Q1(LT(x))+u0(3).*P0(t).*Q2(LT(x)))).^2;
        %I=integral2(u_num,x0,x1,t0,t1);
        I=quad2d(u_num,x0,x1,t0,t1);
    elseif order_s==3
        u_num=@(x,t) (uexact(x,t)-(u0(1).*P0(t).*Q0(x)+u0(2).*P0(t).*Q1(LT(x))+u0(3).*P0(t).*Q2(LT(x))+u0(4).*P0(t).*Q3(LT(x)))).^2;
        %I=integral2(u_num,x0,x1,t0,t1);
        I=quad2d(u_num,x0,x1,t0,t1);
    end

elseif order_t==2
    if order_s==0
       u_num=@(x,t) (uexact(x,t)-(u0.*P0(t).*Q0(x)+u1.*P1(t).*Q0(x))).^2;
       %I=integral2(u_num,x0,x1,t0,t1);
       I=quad2d(u_num,x0,x1,t0,t1);
    elseif order_s==1
       u_num=@(x,t) (uexact(x,t)-(u0(1).*P0(t).*Q0(x)+u0(2).*P0(t).*Q1(LT(x))+u1(1).*P1(t).*Q0(x)+u1(2).*P1(t).*Q1(LT(x)))).^2;
       %I=integral2(u_num,x0,x1,t0,t1);
       I=quad2d(u_num,x0,x1,t0,t1);
    elseif order_s==2
       u_num=@(x,t) (uexact(x,t)-(u0(1).*P0(t).*Q0(x)+u0(2).*P0(t).*Q1(LT(x))+u0(3).*P0(t).*Q2(LT(x))+u1(1).*P1(t).*Q0(x)+u1(2).*P1(t).*Q1(LT(x))+u1(3).*P1(t).*Q2(LT(x)))).^2;
       %I=integral2(u_num,x0,x1,t0,t1);
       I=quad2d(u_num,x0,x1,t0,t1);
    elseif order_s==3
       u_num=@(x,t) (uexact(x,t)-(u0(1).*P0(t).*Q0(x)+u0(2).*P0(t).*Q1(LT(x))+u0(3).*P0(t).*Q2(LT(x))+u0(4).*P0(t).*Q3(LT(x))+u1(1).*P1(t).*Q0(x)+u1(2).*P1(t).*Q1(LT(x))+u1(3).*P1(t).*Q2(LT(x))+u1(4).*P1(t).*Q3(LT(x)))).^2;
       %I=integral2(u_num,x0,x1,t0,t1);
       I=quad2d(u_num,x0,x1,t0,t1);
    end

elseif order_t==3
    if order_s==0
       u_num=@(x,t) (uexact(x,t)-(u0.*P0(t).*Q0(x)+u1.*P1(t).*Q0(x)+u2.*P2(t).*Q0(x))).^2;
       %I=integral2(u_num,x0,x1,t0,t1);
       I=quad2d(u_num,x0,x1,t0,t1);
    elseif order_s==1
       u_num=@(x,t) (uexact(x,t)-(u0(1).*P0(t).*Q0(x)+u0(2).*P0(t).*Q1(LT(x))+u1(1).*P1(t).*Q0(x)+u1(2).*P1(t).*Q1(LT(x))+u2(1).*P2(t).*Q0(x)+u2(2).*P2(t).*Q1(LT(x)))).^2;
       %I=integral2(u_num,x0,x1,t0,t1);
       I=quad2d(u_num,x0,x1,t0,t1);
    elseif order_s==2
       u_num=@(x,t) (uexact(x,t)-(u0(1).*P0(t).*Q0(x)+u0(2).*P0(t).*Q1(LT(x))+u0(3).*P0(t).*Q2(LT(x))+u1(1).*P1(t).*Q0(x)+u1(2).*P1(t).*Q1(LT(x))+u1(3).*P1(t).*Q2(LT(x))+u2(1).*P2(t).*Q0(x)+u2(2).*P2(t).*Q1(LT(x))+u2(3).*P2(t).*Q2(LT(x)))).^2;
       %I=integral2(u_num,x0,x1,t0,t1);
       I=quad2d(u_num,x0,x1,t0,t1);
    elseif order_s==3
       u_num=@(x,t) (uexact(x,t)-(u0(1).*P0(t).*Q0(x)+u0(2).*P0(t).*Q1(LT(x))+u0(3).*P0(t).*Q2(LT(x))+u0(4).*P0(t).*Q3(LT(x))+u1(1).*P1(t).*Q0(x)+u1(2).*P1(t).*Q1(LT(x))+u1(3).*P1(t).*Q2(LT(x))+u1(4).*P1(t).*Q3(LT(x))+u2(1).*P2(t).*Q0(x)+u2(2).*P2(t).*Q1(LT(x))+u2(3).*P2(t).*Q2(LT(x))+u2(4).*P2(t).*Q2(LT(x)))).^2;
       %I=integral2(u_num,x0,x1,t0,t1);
       I=quad2d(u_num,x0,x1,t0,t1);
    end
elseif order_t==4
    if order_s==0
       u_num=@(x,t) (uexact(x,t)-(u0.*P0(t).*Q0(x)+u1.*P1(t).*Q0(x)+u2.*P2(t).*Q0(x)+u3.*P3(t).*Q0(x))).^2;
       %I=integral2(u_num,x0,x1,t0,t1);
       I=quad2d(u_num,x0,x1,t0,t1);
    elseif order_s==1
       u_num=@(x,t) (uexact(x,t)-(u0(1).*P0(t).*Q0(x)+u0(2).*P0(t).*Q1(LT(x))+u1(1).*P1(t).*Q0(x)+u1(2).*P1(t).*Q1(LT(x))+u2(1).*P2(t).*Q0(x)+u2(2).*P2(t).*Q1(LT(x))+u3(1).*P3(t).*Q0(x)+u3(2).*P3(t).*Q1(LT(x)))).^2;
       %I=integral2(u_num,x0,x1,t0,t1);
       I=quad2d(u_num,x0,x1,t0,t1);
    elseif order_s==2
       u_num=@(x,t) (uexact(x,t)-(u0(1).*P0(t).*Q0(x)+u0(2).*P0(t).*Q1(LT(x))+u0(3).*P0(t).*Q2(LT(x))+u1(1).*P1(t).*Q0(x)+u1(2).*P1(t).*Q1(LT(x))+u1(3).*P1(t).*Q2(LT(x))+u2(1).*P2(t).*Q0(x)+u2(2).*P2(t).*Q1(LT(x))+u2(3).*P2(t).*Q2(LT(x))+u3(1).*P3(t).*Q0(x)+u3(2).*P3(t).*Q1(LT(x))+u3(3).*P3(t).*Q2(LT(x)))).^2;
       %I=integral2(u_num,x0,x1,t0,t1);
       I=quad2d(u_num,x0,x1,t0,t1);
    elseif order_s==3
       u_num=@(x,t) (uexact(x,t)-(u0(1).*P0(t).*Q0(x)+u0(2).*P0(t).*Q1(LT(x))+u0(3).*P0(t).*Q2(LT(x))+u0(4).*P0(t).*Q3(LT(x))+u1(1).*P1(t).*Q0(x)+u1(2).*P1(t).*Q1(LT(x))+u1(3).*P1(t).*Q2(LT(x))+u1(4).*P1(t).*Q3(LT(x))+u2(1).*P2(t).*Q0(x)+u2(2).*P2(t).*Q1(LT(x))+u2(3).*P2(t).*Q2(LT(x))+u2(4).*P2(t).*Q2(LT(x))+u3(1).*P3(t).*Q0(x)+u3(2).*P3(t).*Q1(LT(x))+u3(3).*P3(t).*Q2(LT(x))+u3(4).*P3(t).*Q2(LT(x)))).^2;
       %I=integral2(u_num,x0,x1,t0,t1);
       I=quad2d(u_num,x0,x1,t0,t1);
    end
end


end

