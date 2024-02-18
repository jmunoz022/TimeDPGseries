function uin=local_initial(u0,nquad,order,xcoord)
%%%Function that calculates the local initial condition%%%

%%%coordinates of the element%%%
x1=xcoord(1);
x2=xcoord(2);

%%%length of the element%%%
h=x2-x1;        

%%%Quadrature points and weights%%%
[z,w] = quadrature(nquad);

%%%Jacobian%%%
J=abs(h/2);

%%Integration%%%
uin=zeros(order+1,1);
for i=1:nquad
    N=Shape1Q(z(i),order);
    l=(h/2)*z(i)+(x1+x2)/2;
    uin=uin+N*u0(l)*w(i)*J;
end


end