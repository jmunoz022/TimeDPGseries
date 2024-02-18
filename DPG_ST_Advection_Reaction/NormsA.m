%%%Script to study the norms of the matrices%%%
clear all, close all

%%%Data of the problem%%%
[a,b,T,beta,creac,f,g,uin,uexact] = data;

%%Select method in space: method=1 (DPG), method=2 (pDPG)
method=2;
%%Inner product for the Gramm matrix: inner=0 (graph norm) inner=1 (localizable adjoint norm)
inner=1;

%%%Discretization in space%%%
rin=5; rfin=5;
ord_tr=4;
ord_ts=ord_tr+1;
nquad=ord_ts+1;

%%%Discretization in time%%%
order=1;c=0;
k=4;
steps=2^k;
tau=T/steps;
t=0:tau:T;

%%%Efficienty vector%%%
Eff=zeros(rfin-rin+1,1);
i=1;

for r=rin:rfin
    nel=2^r; 
    xsol=linspace(a,b,nel+1);
    
    %%%Select the discretization in space%%%
    [A,MQ,F,A1,A2,A3,A4]=DPG_matrices(beta,creac,xsol,ord_tr,ord_ts,nel,nquad,inner,method);
    dimx=size(A,1);

    %%%L^2 projection of the initial condition%%%
    U0=initial_cond(uin,MQ,xsol,nel,ord_tr,nquad);

    %%%Compute Phi functions%%%
    Phi1=phipade(-tau*A,1);
    Phi2=phipade(-tau*A,2);

    %%%Initialization%%%
    uhat=zeros(nel*(ord_tr+1),steps+1);
    uhat(:,1)=U0;
    uhat1=zeros(nel*(ord_tr+1),steps+1);
    uhat1(:,1)=U0;
    uhat2=zeros(nel*(ord_tr+1),steps+1);
    uhat2(:,1)=U0;
    
    %%%DPG0 method%%%
    for i=1:steps
        %%%source term%%%
        f1=source_term_DPG(f,g,beta,t(i)+c*tau,F,xsol,nel,ord_ts,nquad);

        %%%Compute solution direct method%%%
        uhat(:,i+1)=uhat(:,i)+tau*Phi1*(f1-A*uhat(:,i));

        %%%Compute solution Higham%%%
        W1=f1-A*uhat1(:,i);
        [V1,s1,m1,mv1,mvd1] = expmv(1,[-tau*A,W1;zeros(1,dimx),0],[zeros(dimx,1);1]);
        uhat1(:,i+1)=uhat1(:,i)+tau*V1(1:end-1);
        fprintf('expmv: m = %d, s = %d, prod = %d, prod_Taylor = %d \n', m1, s1, mv1, mvd1)


        %%%Compute solution modified Higham%%%
        U1=A4*uhat2(:,i);
        U=A3*U1;
        W2=f1-A1*uhat2(:,i)+A2*U;
        [V2,s2,m2,mv2,mvd2] = expmv_dpg(1,-tau*A1,-tau*A2,A3,A4,W2,0,[zeros(dimx,1);1]);
        uhat2(:,i+1)=uhat2(:,i)+tau*V2(1:end-1);
        fprintf('expmv_dpg: m = %d, s = %d, prod = %d, prod_Taylor = %d', m2, s2, mv2, mvd2)
        
        norm(uhat(:,i+1)-uhat1(:,i+1))
        norm(uhat(:,i+1)-uhat2(:,i+1))
        pause


    end

end

