function ff=source_term_PG(source,g,beta,t0,F,xsol,nel,order,nquad)

nodes=zeros(nel,2);
for i=1:nel
    nodes(i,1)=i;
    nodes(i,2)=i+1;
end

ff=zeros(order*nel+1,1);
gg=zeros(order*nel+1,1);
gg(1)=1;
for iel=1:nel
    nd=nodes(iel,:);
    xcoord=xsol(nd);
    f=local_source_PG(source,nquad,xcoord,order,t0);
    %%%assembly%%%
    ii=(order*iel-order+1):(order*iel+1);
    ff(ii)=ff(ii)+f;
end  
ff=ff+g(t0)*beta*gg;

%%%homogeneous Dirchlet BC%%%
ff=ff(1:end-1);

ff=F*ff;


end