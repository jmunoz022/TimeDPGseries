function [] = Plots_2D(uhat,udpg0,udpg1,udpg2,udpg3,nel,ord_tr,steps,xsol,t,order)
    %%%Space-time discontinuous problem%%%
    [a,b,T,beta,creac,f,g,uin,uexact] = data;
    
    orderp=ord_tr+1;
    
    %%Plot traces exact solution%%
    % figure
    % hold on
    % for k=1:steps+1
    %     surf([xfine;xfine],repmat(t(k),2,s2),[uexact(xfine,t(k));uexact(xfine,t(k))],'EdgeColor','interp','linewidth',1.5)
    % end
    % colormap(parula)
    % colorbar
    % grid on
    % view(-50,50)
    % xlabel({'$x$'},'interpreter','latex')
    % ylabel({'$t$'},'interpreter','latex')
    % set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',15)
    % hold off
    
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
    
    
    pause

    %%Plot interior%%
    if order==1
        figure
        hold on
        for k=1:steps
            k
            tfinei=linspace(t(k),t(k+1),3);
            s1=size(tfinei,2);
            Uu=reshape(udpg0(:,k),[orderp nel])';
            for i=1:nel
                xfinei=linspace(xsol(i),xsol(i+1),3);
                s2i=size(xfinei,2);
                tfineimat=repmat(tfinei',1,s2i);
                xsolmat=repmat(xfinei,s1,1);
                xl=2*(xfinei-(xsol(i+1)+xsol(i))/2)/(xsol(i+1)-xsol(i));
                usol=zeros(1,size(xl,2));
                for j=1:orderp
                    usol=usol+Uu(i,j)*LegPol(xl,j);
                end
                umati=repmat(usol,s1,1);
                mesh(xsolmat,tfineimat,umati)
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
    elseif order==2 
        figure
        hold on
        for k=1:steps
            tfinei=linspace(t(k),t(k+1),10);
            s1=size(tfinei,2);
            Uu0=reshape(udpg0(:,k),[orderp nel])';
            Uu1=reshape(udpg1(:,k),[orderp nel])';
            for i=1:nel
                xfinei=linspace(xsol(i),xsol(i+1),10);
                s2i=size(xfinei,2);
                tfineimat=repmat(tfinei',1,s2i);
                xsolmat=repmat(xfinei,s1,1);
                xl=2*(xfinei-(xsol(i+1)+xsol(i))/2)/(xsol(i+1)-xsol(i));
                usol0=zeros(1,size(xl,2));
                usol1=zeros(1,size(xl,2));
                for j=1:orderp
                    usol0=usol0+Uu0(i,j)*LegPol(xl,j);
                    usol1=usol1+Uu1(i,j)*LegPol(xl,j);
                end
                umati=repmat(usol0,s1,1)+(usol1'*(tfinei-t(k))/(t(k+1)-t(k)))';
                mesh(xsolmat,tfineimat,umati)
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
    elseif order==3
        figure
        hold on
        for k=1:steps
            tfinei=linspace(t(k),t(k+1),10);
            s1=size(tfinei,2);
            Uu0=reshape(udpg0(:,k),[orderp nel])';
            Uu1=reshape(udpg1(:,k),[orderp nel])';
            Uu2=reshape(udpg2(:,k),[orderp nel])';
            for i=1:nel
                xfinei=linspace(xsol(i),xsol(i+1),10);
                s2i=size(xfinei,2);
                tfineimat=repmat(tfinei',1,s2i);
                xsolmat=repmat(xfinei,s1,1);
                xl=2*(xfinei-(xsol(i+1)+xsol(i))/2)/(xsol(i+1)-xsol(i));
                usol0=zeros(1,size(xl,2));
                usol1=zeros(1,size(xl,2));
                usol2=zeros(1,size(xl,2));
                for j=1:orderp
                    usol0=usol0+Uu0(i,j)*LegPol(xl,j);
                    usol1=usol1+Uu1(i,j)*LegPol(xl,j);
                    usol2=usol2+Uu2(i,j)*LegPol(xl,j);
                end
                umati=repmat(usol0,s1,1)+(usol1'*(tfinei-t(k))/(t(k+1)-t(k)))'+(usol2'*((tfinei-t(k))/(t(k+1)-t(k))).^2)';
                mesh(xsolmat,tfineimat,umati)
            end
        end
        colormap(parula)
        colorbar;
        grid on
        view(-75,25)
        xlabel({'$x$'},'interpreter','latex')
        ylabel({'$t$'},'interpreter','latex')
        set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',15)
        hold off
     elseif order==4
        figure
        hold on
        for k=1:steps
            tfinei=linspace(t(k),t(k+1),10);
            s1=size(tfinei,2);
            Uu0=reshape(udpg0(:,k),[orderp nel])';
            Uu1=reshape(udpg1(:,k),[orderp nel])';
            Uu2=reshape(udpg2(:,k),[orderp nel])';
            Uu3=reshape(udpg3(:,k),[orderp nel])';
            for i=1:nel
                xfinei=linspace(xsol(i),xsol(i+1),10);
                s2i=size(xfinei,2);
                tfineimat=repmat(tfinei',1,s2i);
                xsolmat=repmat(xfinei,s1,1);
                xl=2*(xfinei-(xsol(i+1)+xsol(i))/2)/(xsol(i+1)-xsol(i));
                usol0=zeros(1,size(xl,2));
                usol1=zeros(1,size(xl,2));
                usol2=zeros(1,size(xl,2));
                usol3=zeros(1,size(xl,2));
                for j=1:orderp
                    usol0=usol0+Uu0(i,j)*LegPol(xl,j);
                    usol1=usol1+Uu1(i,j)*LegPol(xl,j);
                    usol2=usol2+Uu2(i,j)*LegPol(xl,j);
                    usol3=usol3+Uu3(i,j)*LegPol(xl,j);
                end
                umati=repmat(usol0,s1,1)+(usol1'*(tfinei-t(k))/(t(k+1)-t(k)))'+(usol2'*((tfinei-t(k))/(t(k+1)-t(k))).^2)'+(usol3'*((tfinei-t(k))/(t(k+1)-t(k))).^3)';
                mesh(xsolmat,tfineimat,umati)
            end
        end
        colormap(parula)
        colorbar
        set(c,'interpreter','latex');
        grid on
        view(-75,25)
        xlabel({'$x$'},'interpreter','latex')
        ylabel({'$t$'},'interpreter','latex')
        set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',15)
        hold off
    end
end




