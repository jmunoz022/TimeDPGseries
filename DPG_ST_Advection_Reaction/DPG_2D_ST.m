%%%DPG method in time for 1D+time advection-reaction equation:%%%
%%%u_t+b*u_x+c*u=f(x,t), u(0,t)=g(t), u(x,0)=u0(x)%%%
%%%after discretizing in space by DPG%%%
%%%U'(t)+A*U(t)=F(t), U(0)=U_0%%%

%%%Data of the problem%%%
[a,b,T,beta,creac,f,g,uin,uexact] = data;

%%Select method in space: method=0 (PG), method=1 (DPG), method=2 (pDPG)
method=2;
%%Inner product for the Gramm matrix: inner=0 (graph norm) inner=1 (localizable adjoint norm)
inner=1;

%%%Discretization in space%%%
nref=2;
ord_tr=1;
ord_ts=ord_tr+1;
nquad=ord_ts+1;

%%%Order and integration points in time%%%
order=1;c=0;
%order=2;c=[0 1];
%order=3;c=[0 1/2 1];
%order=4;c=[0 1/3 2/3 1];

%%%Discretization in time%%%
mref=2;
nsteps=zeros(1,mref+1);

%%%Error%%%
Err_dpg=[];
Err_trace1=[];
Err_trace2=[];
DoF=[];
DoF_trace=[];

for l=1:nref
    nel=2^l; 
    xsol=linspace(a,b,nel+1);
    
    %%%Select the discretization in space%%%
    if method==0
        [A,MQ,F]=PG_matrices(beta,creac,xsol,ord_tr,ord_tr+1,nel,nquad);
    else
        [A,MQ,F]=DPG_matrices(beta,creac,xsol,ord_tr,ord_ts,nel,nquad,inner,method);
    end

    %%%L^2 projection of the initial condition%%%
    U0=initial_cond(uin,MQ,xsol,nel,ord_tr,nquad);
    
    %%%Solve trace variables and interiors%%%
    for k=mref:mref
        steps=2^k;
        nsteps(k+1)=steps;
        tau=T/steps;
        t=0:tau:T;

        if order==1
            %%%Compute Phi functions%%%
            Phi1=phipade(-tau*A,1);
            Phi2=phipade(-tau*A,2);

            %%%Initialization%%%
            uhat=zeros(nel*(ord_tr+1),steps+1);
            uhat(:,1)=U0;
            udpg0=zeros(nel*(ord_tr+1),steps);

            for i=1:steps
                %%%source term%%%
                if method==0
                    f1=source_term_PG(f,g,beta,t(i)+c*tau,F,xsol,nel,ord_ts,nquad);
                else
                    f1=source_term_DPG(f,g,beta,t(i)+c*tau,F,xsol,nel,ord_ts,nquad);
                end
                %%%compute solutions%%%
                uhat(:,i+1)=uhat(:,i)+tau*Phi1*(f1-A*uhat(:,i));
                udpg0(:,i)=Phi1*uhat(:,i)+tau*Phi2*f1;
            end
        
            %%%Plot solution%%%
            %Plots_2D(uhat,udpg0,0,0,0,nel,ord_tr,steps,xsol,t,order);
            Plots_2D_average(uhat,udpg0,0,0,0,nel,ord_tr,steps,xsol,t,order);
        
            %%%%%L2 error of fields%%%%
%             Err=0;
%             for j=1:steps
%                 UU=reshape(udpg0(:,j),[ord_tr+1 nel])';
%                 for i=1:nel
%                     I=Error_fields(uexact,UU(i,:),0,0,0,xsol(i),xsol(i+1),t(j),t(j+1),order,ord_tr);
%                     Err=Err+I;
%                 end
%             end
%             Err_dpg=[Err_dpg sqrt(Err)];
%             DoF=[DoF steps*order*nel*(ord_tr+1)];
% 
%             %%L^2 and Max error of traces%%%
%             Err_trace=zeros(1,steps+1);
%             for j=1:steps+1
%                 Err=0;
%                 UU=reshape(uhat(:,j),[ord_tr+1 nel])';
%                 for i=1:nel
%                     I=Error_trace(uexact,UU(i,:),xsol(i),xsol(i+1),t(j),ord_tr);
%                     Err=Err+I;
%                 end
%                 Err_trace(j)=Err;
%             end
%             Err_trace1=[Err_trace1 sqrt(sum(Err_trace))];
%             Err_trace2=[Err_trace2 max(sqrt(Err_trace))];
%             DoF_trace=[DoF_trace steps*nel*(ord_tr+1)];
        
        elseif order==2
             %%%Compute Phi functions%%%
             Phi1=phipade(-tau*A,1);
             Phi2=phipade(-tau*A,2);
             Phi3=phipade(-tau*A,3);
             Phi4=phipade(-tau*A,4);
             
             %%%Initialization%%%
             uhat=zeros(nel*(ord_tr+1),steps+1);
             uhat(:,1)=U0;
             udpg0=zeros(nel*(ord_tr+1),steps);
             udpg1=zeros(nel*(ord_tr+1),steps);
             
             %%%Compute matrices in time%%%
             matrixinv=kron(invhilb(order),eye(size(A)));
             Atr=[Phi2 Phi1];
             Aint=[Phi3 Phi2;-Phi4+Phi3 -Phi3+Phi2];
             L=LgrC(c,order);
             B=kron(L,eye(size(A))); 

             for i=1:steps
                %%%Source term%%% 
                if method==0
                    f1=source_term_PG(f,g,beta,t(i)+c(1)*tau,F,xsol,nel,ord_ts,nquad);
                    f2=source_term_PG(f,g,beta,t(i)+c(2)*tau,F,xsol,nel,ord_ts,nquad);
                else
                    f1=source_term_DPG(f,g,beta,t(i)+c(1)*tau,F,xsol,nel,ord_ts,nquad);
                    f2=source_term_DPG(f,g,beta,t(i)+c(2)*tau,F,xsol,nel,ord_ts,nquad);
                end
                %%%Compute solution%%%
                uhat(:,i+1)=uhat(:,i)-tau*Phi1*A*uhat(:,i)+tau*Atr*B*[f1;f2];
                rhs=[Phi1;Phi1-Phi2]*uhat(:,i)+tau*Aint*B*[f1;f2];
                usol=matrixinv*rhs;
                udpg0(:,i)=usol(1:nel*(ord_tr+1),:);
                udpg1(:,i)=usol((nel*(ord_tr+1)+1):end,:);
             end
             
             %%%%Plot solution%%%
             %Plots_2D(uhat,udpg0,udpg1,0,0,nel,ord_tr,steps,xsol,t,order);
             Plots_2D_average(uhat,udpg0,udpg1,0,0,nel,ord_tr,steps,xsol,t,order);
             
             %%%L2 error of fields%%%
             Err=0;
             for j=1:steps
                 UU0=reshape(udpg0(:,j),[ord_tr+1 nel])';
                 UU1=reshape(udpg1(:,j),[ord_tr+1 nel])';
                 for i=1:nel
                     I=Error_fields(uexact,UU0(i,:),UU1(i,:),0,0,xsol(i),xsol(i+1),t(j),t(j+1),order,ord_tr);
                     Err=Err+I;
                 end
             end
             Err_dpg=[Err_dpg sqrt(Err)];
             DoF=[DoF steps*order*nel*(ord_tr+1)];

             %%%L^2 and Max error of traces%%%
             Err_trace=zeros(1,steps+1);
             for j=1:steps+1
                 Err=0;
                 UU=reshape(uhat(:,j),[ord_tr+1 nel])';
                 for i=1:nel
                     I=Error_trace(uexact,UU(i,:),xsol(i),xsol(i+1),t(j),ord_tr);
                     Err=Err+I;
                 end
                 Err_trace(j)=Err;
             end
             Err_trace1=[Err_trace1 sqrt(sum(Err_trace))];
             Err_trace2=[Err_trace2 max(sqrt(Err_trace))];
             DoF_trace=[DoF_trace steps*nel*(ord_tr+1)];
        
        elseif order==3
             %%%Compute Phi functions%%%
             Phi1=phipade(-tau*A,1);
             Phi2=phipade(-tau*A,2);
             Phi3=phipade(-tau*A,3);
             Phi4=phipade(-tau*A,4);
             Phi5=phipade(-tau*A,5);
             Phi6=phipade(-tau*A,6);
             
             %%%Initialization%%%
             uhat=zeros(nel*(ord_tr+1),steps+1);
             uhat(:,1)=U0;
             udpg0=zeros(nel*(ord_tr+1),steps);
             udpg1=zeros(nel*(ord_tr+1),steps);
             udpg2=zeros(nel*(ord_tr+1),steps);
             
             %%%Compute matrices in time%%%
             matrixinv=kron(invhilb(order),eye(size(A)));
             Atr=[2*Phi3 Phi2 Phi1];
             Aint=[2*Phi4 Phi3 Phi2;-2*Phi5+2*Phi4 -Phi4+Phi3 -Phi3+Phi2; 4*Phi6-4*Phi5+2*Phi4 2*Phi5-2*Phi4+Phi3 2*Phi4-2*Phi3+Phi2];
             L=LgrC(c,order);
             B=kron(L,eye(size(A))); 

             for i=1:steps
                %%%Source term%%%
                if method==0
                    f1=source_term_PG(f,g,beta,t(i)+c(1)*tau,F,xsol,nel,ord_ts,nquad);
                    f2=source_term_PG(f,g,beta,t(i)+c(2)*tau,F,xsol,nel,ord_ts,nquad);
                    f3=source_term_PG(f,g,beta,t(i)+c(3)*tau,F,xsol,nel,ord_ts,nquad);
                else
                    f1=source_term_DPG(f,g,beta,t(i)+c(1)*tau,F,xsol,nel,ord_ts,nquad);
                    f2=source_term_DPG(f,g,beta,t(i)+c(2)*tau,F,xsol,nel,ord_ts,nquad);
                    f3=source_term_DPG(f,g,beta,t(i)+c(3)*tau,F,xsol,nel,ord_ts,nquad);
                end
                %%%Compute solution%%%
                uhat(:,i+1)=uhat(:,i)-tau*Phi1*A*uhat(:,i)+tau*Atr*B*[f1;f2;f3];
                rhs=[Phi1;Phi1-Phi2;Phi1-2*Phi2+2*Phi3]*uhat(:,i)+tau*Aint*B*[f1;f2;f3];
                usol=matrixinv*rhs;
                udpg0(:,i)=usol(1:nel*(ord_tr+1),:);
                udpg1(:,i)=usol((nel*(ord_tr+1)+1):2*nel*(ord_tr+1));
                udpg2(:,i)=usol((2*nel*(ord_tr+1)+1):end);
             end
            
             %%%%Plot solution%%%
             %Plots_2D(uhat,udpg0,udpg1,udpg2,0,nel,ord_tr,steps,xsol,t,order);
             %Plots_2D_average(uhat,udpg0,udpg1,udpg2,0,nel,ord_tr,steps,xsol,t,order);
             Plots_1D
             
             %%%L2 error of fields%%%
             Err=0;
             for j=1:steps
                 UU0=reshape(udpg0(:,j),[ord_tr+1 nel])';
                 UU1=reshape(udpg1(:,j),[ord_tr+1 nel])';
                 UU2=reshape(udpg2(:,j),[ord_tr+1 nel])';
                 for i=1:nel
                     I=Error_fields(uexact,UU0(i,:),UU1(i,:),UU2(i,:),0,xsol(i),xsol(i+1),t(j),t(j+1),order,ord_tr);
                     Err=Err+I;
                 end
             end
             Err_dpg=[Err_dpg sqrt(Err)];
             DoF=[DoF steps*order*nel*(ord_tr+1)];

             %%%L^2 and Max error of traces%%%
             Err_trace=zeros(1,steps+1);
             for j=1:steps+1
                 Err=0;
                 UU=reshape(uhat(:,j),[ord_tr+1 nel])';
                 for i=1:nel
                     I=Error_trace(uexact,UU(i,:),xsol(i),xsol(i+1),t(j),ord_tr);
                     Err=Err+I;
                 end
                 Err_trace(j)=Err;
             end
             Err_trace1=[Err_trace1 sqrt(sum(Err_trace))];
             Err_trace2=[Err_trace2 max(sqrt(Err_trace))];
             DoF_trace=[DoF_trace steps*nel*(ord_tr+1)];   

        elseif order==4
             %%%Compute Phi functions%%%
             Phi1=phipade(-tau*A,1);
             Phi2=phipade(-tau*A,2);
             Phi3=phipade(-tau*A,3);
             Phi4=phipade(-tau*A,4);
             Phi5=phipade(-tau*A,5);
             Phi6=phipade(-tau*A,6);
             Phi7=phipade(-tau*A,7);
             Phi8=phipade(-tau*A,8);
                         
             %%%Initialization%%%
             uhat=zeros(nel*(ord_tr+1),steps+1);
             uhat(:,1)=U0;
             udpg0=zeros(nel*(ord_tr+1),steps);
             udpg1=zeros(nel*(ord_tr+1),steps);
             udpg2=zeros(nel*(ord_tr+1),steps);
             udpg3=zeros(nel*(ord_tr+1),steps);
             
             %%%Compute matrices in time%%%
             matrixinv=kron(invhilb(order),eye(size(A)));
             Atr=[6*Phi4 2*Phi3 Phi2 Phi1];
             Aint=[6*Phi5 2*Phi4 Phi3 Phi2;-6*Phi6+6*Phi5 -2*Phi5+2*Phi4 -Phi4+Phi3 -Phi3+Phi2; 12*Phi7-12*Phi6+6*Phi5 4*Phi6-4*Phi5+2*Phi4 2*Phi5-2*Phi4+Phi3 2*Phi4-2*Phi3+Phi2; -36*Phi8+36*Phi7-18*Phi6+6*Phi5 -12*Phi7+12*Phi6-6*Phi5+2*Phi4 -6*Phi6+6*Phi5-3*Phi4+Phi3 -6*Phi5+6*Phi4-3*Phi3+Phi2];
             L=LgrC(c,order);
             B=kron(L,eye(size(A)));  
             
             for i=1:steps
                %%%Source term%%%
                if method==0
                    f1=source_term_PG(f,g,beta,t(i)+c(1)*tau,F,xsol,nel,ord_ts,nquad);
                    f2=source_term_PG(f,g,beta,t(i)+c(2)*tau,F,xsol,nel,ord_ts,nquad);
                    f3=source_term_PG(f,g,beta,t(i)+c(3)*tau,F,xsol,nel,ord_ts,nquad);
                    f4=source_term_PG(f,g,beta,t(i)+c(4)*tau,F,xsol,nel,ord_ts,nquad);
                else
                    f1=source_term_DPG(f,g,beta,t(i)+c(1)*tau,F,xsol,nel,ord_ts,nquad);
                    f2=source_term_DPG(f,g,beta,t(i)+c(2)*tau,F,xsol,nel,ord_ts,nquad);
                    f3=source_term_DPG(f,g,beta,t(i)+c(3)*tau,F,xsol,nel,ord_ts,nquad);
                    f4=source_term_DPG(f,g,beta,t(i)+c(4)*tau,F,xsol,nel,ord_ts,nquad);
                end
                %%%Compute solution%%%
                uhat(:,i+1)=uhat(:,i)-tau*Phi1*A*uhat(:,i)+tau*Atr*B*[f1;f2;f3;f4];
                rhs=[Phi1;Phi1-Phi2;Phi1-2*Phi2+2*Phi3; Phi1-3*Phi2+6*Phi3-6*Phi4]*uhat(:,i)+tau*Aint*B*[f1;f2;f3;f4];
                usol=matrixinv*rhs;
                udpg0(:,i)=usol(1:nel*(ord_tr+1),:);
                udpg1(:,i)=usol((nel*(ord_tr+1)+1):2*nel*(ord_tr+1));
                udpg2(:,i)=usol((2*nel*(ord_tr+1)+1):3*nel*(ord_tr+1));
                udpg3(:,i)=usol((3*nel*(ord_tr+1)+1):end);
             end
            
             %%%%Plot solution%%%
             %Plots_2D(uhat,udpg0,udpg1,udpg2,udpg3,nel,ord_tr,steps,xsol,t,order);
             
             %%%L^2 error of fields%%%
             Err=0;
             for j=1:steps
                 UU0=reshape(udpg0(:,j),[ord_tr+1 nel])';
                 UU1=reshape(udpg1(:,j),[ord_tr+1 nel])';
                 UU2=reshape(udpg2(:,j),[ord_tr+1 nel])';
                 UU3=reshape(udpg3(:,j),[ord_tr+1 nel])';
                 for i=1:nel
                     I=Error_fields(uexact,UU0(i,:),UU1(i,:),UU2(i,:),UU3(i,:),xsol(i),xsol(i+1),t(j),t(j+1),order,ord_tr);
                     Err=Err+I;
                 end
             end
             Err_dpg=[Err_dpg sqrt(Err)];
             DoF=[DoF steps*order*nel*(ord_tr+1)];

             %%%L^2 and Max error of traces%%%
             Err_trace=zeros(1,steps+1);
             for j=1:steps+1
                 Err=0;
                 UU=reshape(uhat(:,j),[ord_tr+1 nel])';
                 for i=1:nel
                     I=Error_trace(uexact,UU(i,:),xsol(i),xsol(i+1),t(j),ord_tr);
                     Err=Err+I;
                 end
                 Err_trace(j)=Err;
             end
             Err_trace1=[Err_trace1 sqrt(sum(Err_trace))];
             Err_trace2=[Err_trace2 max(sqrt(Err_trace))];
             DoF_trace=[DoF_trace steps*nel*(ord_tr+1)];        
        end
    end
end

%Plot convergence curves fields
figure
hold on
plot(DoF,Err_dpg,'k-o')
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel({'$DoF$'},'interpreter','latex')
ylabel({'$||u-u_h||_U$'},'interpreter','latex')
set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',15)
hold off

%Plot convergence curves traces
figure
hold on
plot(DoF_trace,Err_trace1,'b-*',DoF_trace,Err_trace2,'k-o')
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel({'$DoF$'},'interpreter','latex')
ylabel({'$||u-u_h||_U$'},'interpreter','latex')
set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',15)
hold off

% %Save error in a .txt file
dof=DoF';
error=Err_dpg';
Tab=table(dof,error)
writetable(Tab,'ConvergenceFieldsTime2.txt','Delimiter',' ');
% 
% dof=DoF_trace';
% error=Err_trace1';
% Tab=table(dof,error)
% writetable(Tab,'ConvergenceTraceL2Space2.txt','Delimiter',' ');
% 
% dof=DoF_trace';
% error=Err_trace2';
% Tab=table(dof,error)
% writetable(Tab,'ConvergenceTraceMaxSpace2.txt','Delimiter',' ');
