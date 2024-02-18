function Uin=initial_cond(u0,MQ,xsol,nel,order,nquad)

nodes=zeros(nel,2);
for i=1:nel
    nodes(i,1)=i;
    nodes(i,2)=i+1;
end

uuin=zeros(nel*(order+1),1);
for iel=1:nel
    nd=nodes(iel,:);
    xcoord=xsol(nd);
    uin=local_initial(u0,nquad,order,xcoord);
    %%%assembly%%%
    ii=((order+1)*iel-(order+1)+1):((order+1)*iel);
    uuin(ii)=uuin(ii)+uin;
end  

%%%homogeneous Dirchlet BC%%%
Uin=MQ\uuin;

end