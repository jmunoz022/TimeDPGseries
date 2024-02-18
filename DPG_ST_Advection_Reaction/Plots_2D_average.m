function [] = Plots_2D_average(uhat,udpg0,udpg1,udpg2,udpg3,nel,ord_tr,steps,xsol,t,order)
    %%%Space-time discontinuous problem%%%
    [a,b,T,beta,creac,f,g,uin,uexact] = data;
    
    orderp=ord_tr+1;
    
    tau=t(2)-t(1);

    %%Plot traces%%
    figure
    hold on
    for k=1:steps+1
        k
        Uu=reshape(uhat(:,k),[orderp nel])';
        for i=1:nel
            xfinei=linspace(xsol(i),xsol(i+1),3);
            s2i=size(xfinei,2);
            xl=2*(xfinei-(xsol(i+1)+xsol(i))/2)/(xsol(i+1)-xsol(i));
            usol=zeros(1,size(xl,2));
            for j=1:orderp
                usol=usol+Uu(i,j)*LegPol(xl,j);
            end
            surf([xfinei;xfinei],repmat(t(k),2,s2i),[usol;usol],'EdgeColor','interp','linewidth',1.5)
        end
    end
    colormap(parula)
    colorbar
    grid on
    view(-75,25)
    xlabel({'$x$'},'interpreter','latex')
    ylabel({'$t$'},'interpreter','latex')
    set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',15)
    hold off

    %%Plot averages%%
    figure
    hold on
    for k=1:steps
        k
        if order==1
            Uu=reshape(udpg0(:,k),[orderp nel])';
            for i=1:nel
                xfinei=linspace(xsol(i),xsol(i+1),3);
                s2i=size(xfinei,2);
                xl=2*(xfinei-(xsol(i+1)+xsol(i))/2)/(xsol(i+1)-xsol(i));
                usol=zeros(1,size(xl,2));
                for j=1:orderp
                    usol=usol+Uu(i,j)*LegPol(xl,j);
                end
                surf([xfinei;xfinei],repmat(t(k)+tau/2,2,s2i),[usol;usol],'EdgeColor','interp','linewidth',1.5)
            end
        elseif order==2
            Uu=reshape(udpg0(:,k)+udpg1(:,k)/2,[orderp nel])';
            for i=1:nel
                xfinei=linspace(xsol(i),xsol(i+1),3);
                s2i=size(xfinei,2);
                xl=2*(xfinei-(xsol(i+1)+xsol(i))/2)/(xsol(i+1)-xsol(i));
                usol=zeros(1,size(xl,2));
                for j=1:orderp
                    usol=usol+Uu(i,j)*LegPol(xl,j);
                end
                surf([xfinei;xfinei],repmat(t(k)+tau/2,2,s2i),[usol;usol],'EdgeColor','interp','linewidth',1.5)
            end
        elseif order==3
            Uu=reshape(udpg0(:,k)+udpg1(:,k)/2+udpg2(:,k)/3,[orderp nel])';
            for i=1:nel
                xfinei=linspace(xsol(i),xsol(i+1),3);
                s2i=size(xfinei,2);
                xl=2*(xfinei-(xsol(i+1)+xsol(i))/2)/(xsol(i+1)-xsol(i));
                usol=zeros(1,size(xl,2));
                for j=1:orderp
                    usol=usol+Uu(i,j)*LegPol(xl,j);
                end
                surf([xfinei;xfinei],repmat(t(k)+tau/2,2,s2i),[usol;usol],'EdgeColor','interp','linewidth',1.5)
            end
        elseif order==4
            Uu=reshape(udpg0(:,k)+udpg1(:,k)/2+udpg2(:,k)/3+udpg3(:,k)/4,[orderp nel])';
            for i=1:nel
                xfinei=linspace(xsol(i),xsol(i+1),3);
                s2i=size(xfinei,2);
                xl=2*(xfinei-(xsol(i+1)+xsol(i))/2)/(xsol(i+1)-xsol(i));
                usol=zeros(1,size(xl,2));
                for j=1:orderp
                    usol=usol+Uu(i,j)*LegPol(xl,j);
                end
                surf([xfinei;xfinei],repmat(t(k)+tau/2,2,s2i),[usol;usol],'EdgeColor','interp','linewidth',1.5)
            end
        end
    end
    colormap(parula)
    colorbar
    grid on
    view(-75,25)
    xlabel({'$x$'},'interpreter','latex')
    ylabel({'$t$'},'interpreter','latex')
    set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',15)
    hold off
    
    
end




