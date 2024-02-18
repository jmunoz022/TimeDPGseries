function [val,dval] = Shape1H(x,order)
n=size(x,2);
P=zeros(order+1,n);
dP=zeros(order+1,n);
if order==1
    P(1,:)=(ones(1,n)-x)/2;
    P(end,:)=(x+ones(1,n))/2;
    dP(1,:)=-ones(1,n)/2;
    dP(end,:)=ones(1,n)/2;
    val=P;
    dval=dP;
else 
    P(1,:)=(ones(1,n)-x)/2;
    P(end,:)=(x+ones(1,n))/2;
    dP(1,:)=-ones(1,n)/2;
    dP(end,:)=ones(1,n)/2;
    for i=2:order
        [Pol0,dPol0]=Shape1Q(x,i-2);
        [Pol1,dPol1]=Shape1Q(x,i-1);
        P(i,:)=(i-1)*(Pol0(end,:)-x.*Pol1(end,:));
        dP(i,:)=(i-1)*(dPol0(end,:)-Pol1(end,:)-x.*dPol1(end,:));
    end
    %val=P(order,:);
    %dval=dP(order,:);
    val=P;
    dval=dP;
end
end