function [A,MQ,F]=PG_matrices(beta,creac,xsol,order_trial,order_test,nel,nquad)

nodes=zeros(nel,2);
for i=1:nel
    nodes(i,1)=i;
    nodes(i,2)=i+1;
end

%%%mass and stiffness matrix%%%
kk=zeros((order_trial+1)*nel,order_test*nel+1); 
mm=zeros((order_trial+1)*nel,order_test*nel+1);
mmc=zeros((order_trial+1)*nel,order_test*nel+1);
mmq=zeros((order_trial+1)*nel,(order_trial+1)*nel); 

for iel=1:nel
    nd=nodes(iel,:);
    xcoord=xsol(nd);

    m=local_massH(xcoord,order_trial,order_test,1,nquad);
    mc=local_massH(xcoord,order_trial,order_test,creac,nquad);
    k=local_stiffnessH(xcoord,order_trial,order_test,beta,nquad);
    mq=local_massq(xcoord,order_trial,nquad);


    %%%assembly%%%
    ii=((order_trial+1)*iel-(order_trial+1)+1):((order_trial+1)*iel);
    jj=(order_test*iel-order_test+1):(order_test*iel+1);
    mm(ii,jj)=mm(ii,jj)+m;
    mmc(ii,jj)=mmc(ii,jj)+mc;
    kk(ii,jj)=kk(ii,jj)+k;
    mmq(ii,ii)=mmq(ii,ii)+mq;

end
MQ=mmq';
M=mm';
K=-kk'+mmc';

%%%Petrov-Galerkin%%%
M=M(1:end-1,:);
K=K(1:end-1,:);
A=M\K;
F=inv(M);


%%%Petrov-Galerkin with one trace%%%
% gg=zeros(order_test*nel+1,order_test*nel+1);
% for i=1:nel
%     g=local_gramm(xcoord,order_test,nquad,inner);
%     gg(jj,jj)=gg(jj,jj)+g;
% end
% G0=gg';
% Ghat=zeros(size(G0));
% Ghat(end,end)=1;
% G=G0+Ghat;
% 
% L=zeros(size(G,1),1);
% L(end)=1;
% 
% A0=inv(K'*inv(G)*M)*(K'*inv(G)*K);
% A1=inv(K'*inv(G)*M-K'*inv(G)*L*inv(L'*inv(G)*L)*L'*inv(G)*M)*(K'*inv(G)*K-K'*inv(G)*L*inv(L'*inv(G)*L)*L'*inv(G)*K);

end




