function [A,MQ,F,A1,A2,A3,A4]=DPG_matrices(beta,creac,xsol,order_trial,order_test,nel,nquad,inner,method)

nodes=zeros(nel,2);
for i=1:nel
    nodes(i,1)=i;
    nodes(i,2)=i+1;
end

%%%mass and stiffness matrix%%%
kk=zeros((order_trial+1)*nel,(order_test+1)*nel);
mm=zeros((order_trial+1)*nel,(order_test+1)*nel);
mmc=zeros((order_trial+1)*nel,(order_test+1)*nel);
mmq=zeros((order_trial+1)*nel,(order_trial+1)*nel); 
gg=zeros((order_test+1)*nel,(order_test+1)*nel);
ll=zeros((order_test+1)*nel,nel+1);

for iel=1:nel
    nd=nodes(iel,:);
    xcoord=xsol(nd);

    m=local_mass(xcoord,order_trial,order_test,1,nquad);
    mc=local_mass(xcoord,order_trial,order_test,creac,nquad);
    k=local_stiffness(xcoord,order_trial,order_test,beta,nquad);
    mq=local_massq(xcoord,order_trial,nquad);
    g=local_gramm(xcoord,order_test,nquad,inner);
    l=local_trace(order_test);

    %%%assembly%%%
    ii=((order_trial+1)*iel-(order_trial+1)+1):((order_trial+1)*iel);
    jj=((order_test+1)*iel-(order_test+1)+1):((order_test+1)*iel);
    qq=iel:iel+1;
    mm(ii,jj)=mm(ii,jj)+m;
    mmc(ii,jj)=mmc(ii,jj)+mc;
    kk(ii,jj)=kk(ii,jj)+k;
    mmq(ii,ii)=mmq(ii,ii)+mq;
    gg(jj,jj)=gg(jj,jj)+g;
    ll(jj,qq)=ll(jj,qq)+l;
end
MQ=mmq';
M=mm';
K=-kk'+mmc';
L=ll(:,2:end);
G=gg';

invG=inv(G);

if method==1
    S1=invG-invG*L*inv(L'*invG*L)*L'*invG;
    F=inv(K'*S1*M)*(K'*S1);
    A=F*K;
    
elseif method==2
    S2=invG-invG*M*inv(K'*invG*M)*K'*invG;
    S3=invG-invG*L*inv(L'*S2*L)*L'*S2;
    F=inv(K'*invG*M)*K'*S3;
    A=F*K;

    A1=inv(K'*invG*M)*K'*invG*K;
    A2=inv(K'*invG*M)*K'*invG*L;
    A3=inv(L'*S2*L);
    A4=L'*S2*K;

end


end