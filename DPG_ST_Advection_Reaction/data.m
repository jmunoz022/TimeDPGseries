function [a,b,T,beta,creac,f,g,uin,uexact] = data
%Function that selects the problem
P=4; %parameter to change the problem
%Space-time intervals
T=1;
a=0;
b=1;
if P==1 %Discontinuous transport
    beta=1;
    creac=0;
    uin=@(x) -(x.^2).*(x>=0.5)+(x.^2).*(x<0.5);
    f=@(x,t) 0.*x*t;
    g=@(t) 0.*t;
    ub=@(t) 0.*t;
    uexact=@(x,t) uin(x-t).*(x>=t)+ub(t-x).*(x<t);
elseif P==2 %Traveling initial condition
    beta=1;
    creac=0;
    %uin=@(x) exp(-300*(2*x-0.3).^2).*(abs(2*x-0.3)<=0.25)+1.*(abs(2*x-0.9)<=0.2)+sqrt((1-((2*x-1.6)/(0.2)).^2)).*(abs(2*x-1.6)<=0.2);
    uin=@(x) ((x-0.25)/0.25).*(x>=0.25).*(x<=0.5)+((0.75-x)/0.25).*(x>=0.5).*(x<=0.75);
    f=@(x,t) 0.*x.*t;
    g=@(t) uin(-t);
    uexact=@(x,t) uin(x-t);
elseif P==3 %Smooth convection
    beta=1;
    creac=0;
    f=@(x,t) creac*(1+x.^2+t.^2)+2*t+beta*(2*x);
    g=@(t) 1+t.^2;
    uin=@(x) 1+x.^2;
    uexact=@(x,t) 1+x.^2+t.^2;
elseif P==4 %Smooth convection-reaction
    beta=1;
    creac=0;
    f=@(x,t) beta*pi*cos(pi*x).*cos(pi*t)-pi*sin(pi*x).*sin(pi*t)+creac*sin(pi*x).*cos(pi*t);
    g=@(t) 0;
    uin=@(x) sin(pi*x);
    uexact=@(x,t) sin(pi*x).*cos(pi*t);
elseif P==5 %Smooth convection-reaction
    beta=1;
    creac=1;
    f=@(x,t) (sin(pi*(-t+x)).^2-(-1+t).*sin(pi*(-t+x)).^4)./(-1+t.*sin(pi*(-t+x)).^2).^2;
    g=@(t)  (sin(pi*(-t)).^2)./(1-t*sin(pi*(-t)).^2);
    uin=@(x) (sin(pi*(x)).^2);
    uexact=@(x,t) (sin(pi*(x-t)).^2)./(1-t.*sin(pi*(x-t)).^2);
end

end