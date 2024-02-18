%%%Space-time regular problem%%%
xfine=linspace(a,b,nel*100);
tc=0.25;
k=find(t==tc);
tm=tc+tau*0.5;

ordp=ord_tr+1;

%Plot trace in time 
figure
hold on
U=uhat(:,k);
U=reshape(U,[ordp nel])';
plot(xfine,uexact(xfine,t(k)),'red','LineWidth',1.5)
for i=1:nel
   xfinei=linspace(xsol(i),xsol(i+1),100);
   xl=2*(xfinei-(xsol(i+1)+xsol(i))/2)/(xsol(i+1)-xsol(i));
   usol=zeros(1,size(xl,2));
   for j=1:ordp
       usol=usol+U(i,j)*LegPol(xl,j);
   end
   plot(xfinei,usol,'blue','LineWidth',1.5);
end
xlabel({'$x$'},'interpreter','latex')
set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',20)
title({['$t=\,$',num2str(t(k))]},'interpreter','latex','fontsize',20)
legend({'$u$','$\hat{u}_h$'},'Location','northwest','interpreter','latex')
box on
hold off


% %Plot element interior in time 
% if order==1
%     figure
%     hold on
%     U=udpg0(:,k);
%     U=reshape(U,[ordp nel])';
%     plot(xfine,uexact(xfine,tm),'red','LineWidth',1.5)
%     for i=1:nel
%        xfinei=linspace(xsol(i),xsol(i+1),100);
%        xl=2*(xfinei-(xsol(i+1)+xsol(i))/2)/(xsol(i+1)-xsol(i));
%        usol=zeros(1,size(xl,2));
%        for j=1:ordp
%            usol=usol+U(i,j)*LegPol(xl,j);
%        end
%        plot(xfinei,usol,'blue','LineWidth',1.5);
%     end
%     xlabel({'$x$'},'interpreter','latex')
%     set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',20)
%     title({['$t=\,$',num2str(tm)]},'interpreter','latex','fontsize',20)
%     legend({'$u$','$u_h$'},'Location','northwest','interpreter','latex')
%     box on
%     hold off
% elseif order==2
%     figure
%     hold on
%     U=udpg0(:,k)+udpg1(:,k).*(tm-t(k))/(t(k+1)-t(k));
%     U=reshape(U,[ordp nel])';
%     plot(xfine,uexact(xfine,tm),'red','LineWidth',1.5)
%     for i=1:nel
%        xfinei=linspace(xsol(i),xsol(i+1),100);
%        xl=2*(xfinei-(xsol(i+1)+xsol(i))/2)/(xsol(i+1)-xsol(i));
%        usol=zeros(1,size(xl,2));
%        for j=1:ordp
%            usol=usol+U(i,j)*LegPol(xl,j);
%        end
%        plot(xfinei,usol,'blue','LineWidth',1.5);
%     end
%     xlabel({'$x$'},'interpreter','latex')
%     set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',20)
%     title({['$t=\,$',num2str(tm)]},'interpreter','latex','fontsize',20)
%     %legend({'$u$','$u_h$'},'Location','northwest','interpreter','latex')
%     box on
%     hold off
% elseif order==3
%     figure
%     hold on
%     U=udpg0(:,k)+udpg1(:,k).*(tm-t(k))/(t(k+1)-t(k))+udpg2(:,k).*((tm-t(k))/(t(k+1)-t(k))).^2;
%     U=reshape(U,[ordp nel])';
%     plot(xfine,uexact(xfine,tm),'red','LineWidth',1.5)
%     for i=1:nel
%        xfinei=linspace(xsol(i),xsol(i+1),100);
%        xl=2*(xfinei-(xsol(i+1)+xsol(i))/2)/(xsol(i+1)-xsol(i));
%        usol=zeros(1,size(xl,2));
%        for j=1:ordp
%            usol=usol+U(i,j)*LegPol(xl,j);
%        end
%        plot(xfinei,usol,'blue','LineWidth',1.5);
%     end
%     xlabel({'$x$'},'interpreter','latex')
%     set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',20)
%     title({['$t=\,$',num2str(tm)]},'interpreter','latex','fontsize',20)
%     %legend({'$u$','$u_h$'},'Location','northwest','interpreter','latex')
%     box on
%     hold off
% elseif order==4
%     figure
%     hold on
%     U=udpg0(:,k)+udpg1(:,k).*(tm-t(k))/(t(k+1)-t(k))+udpg2(:,k).*((tm-t(k))/(t(k+1)-t(k))).^2+udpg3(:,k).*((tm-t(k))/(t(k+1)-t(k))).^3;
%     U=reshape(U,[ordp nel])';
%     plot(xfine,uexact(xfine,tm),'red','LineWidth',1.5)
%     for i=1:nel
%        xfinei=linspace(xsol(i),xsol(i+1),100);
%        xl=2*(xfinei-(xsol(i+1)+xsol(i))/2)/(xsol(i+1)-xsol(i));
%        usol=zeros(1,size(xl,2));
%        for j=1:ordp
%            usol=usol+U(i,j)*LegPol(xl,j);
%        end
%        plot(xfinei,usol,'blue','LineWidth',1.5);
%     end
%     xlabel({'$x$'},'interpreter','latex')
%     set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',20)
%     title({['$t=\,$',num2str(tm)]},'interpreter','latex','fontsize',20)
%     %legend({'$u$','$u_h$'},'Location','northwest','interpreter','latex')
%     box on
%     hold off
% end

%savefig('Dis_p2_elem128_tau128_t025.fig')
    
