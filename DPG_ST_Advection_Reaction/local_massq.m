function m=local_massq(xcoord,order,nquad)
%%%Function that stores the local mass matrix for each element%%% 

%xcoord=x coordinates of the element
%nquad=number of quadrature points

x1=xcoord(1);
x2=xcoord(2);
h=x2-x1;

%%%Quadrature points and weights%%%
[z,w] = quadrature(nquad);

%%%Jacobian%%%
J=h/2;

%%%Integration%%%
m=zeros(order+1,order+1);
for i=1:nquad
    N=Shape1Q(z(i),order);
    m=m+N*N'*w(i)*abs(J);
end

end

