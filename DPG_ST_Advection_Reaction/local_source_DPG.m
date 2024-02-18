function fn=local_source_DPG(source,nquad,xcoord,order,t0)
%%%Function that calculates the local source vector at a fixed time tn%%%

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
fn=zeros(order+1,1);
for i=1:nquad
    N=Shape1Q(z(i),order);
    l=(h/2)*z(i)+(x1+x2)/2;
    fn=fn+N*source(l,t0)*w(i)*J;
end


end