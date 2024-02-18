function k=local_stiffness(xcoord,order_trial,order_test,beta,nquad)
%%%Function that stores the local stiffness matrix for each element%%% 

%xcoord=x coordinates of the element
%nquad=number of quadrature points
%diff=diffusion coefficient

x1=xcoord(1);
x2=xcoord(2);
h=x2-x1;

%%%Quadrature points and weights%%%
[z,w] = quadrature(nquad);

%%%Jacobian%%%
J=h/2;

%%%Integration%%%
k=zeros(order_trial+1,order_test+1);
for i=1:nquad
    [N,dN]=Shape1Q(z(i),order_trial);
    [M,dM]=Shape1Q(z(i),order_test);
    k=k+N*dM'*beta*w(i)*abs(J^-1)*abs(J);
end

end

