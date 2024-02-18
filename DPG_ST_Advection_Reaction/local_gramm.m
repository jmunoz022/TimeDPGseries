function g=local_gramm(xcoord,order,nquad,inner)
%%%Function that stores the local gramm matrix for each element%%% 

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
g=zeros(order+1,order+1);

if inner==0
    for i=1:nquad
        [M,dM]=Shape1Q(z(i),order);
        g=g+dM*dM'*w(i)*abs(J^-1)+M*M'*w(i)*abs(J);
    end
elseif inner==1
    for i=1:nquad
        [M,dM]=Shape1Q(z(i),order);
        g=g+dM*dM'*w(i)*abs(J^-1);
    end
    [M,dM]=Shape1Q(1,order);
    g=g+M*M';
end

end

