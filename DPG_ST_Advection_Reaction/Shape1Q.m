function [val,dval] = Shape1Q(x,order)
n=size(x,2);
P=zeros(order+1,n);
dP=zeros(order+1,n);
if order==0
    P(1,:)=ones(1,n);
    dP(1,:)=zeros(1,n);
elseif order==1
    P(1,:)=ones(1,n);
    dP(1,:)=zeros(1,n);
    P(2,:)=x;
    dP(2,:)=ones(1,n);
else 
    P(1,:)=ones(1,n);
    dP(1,:)=zeros(1,n);
    P(2,:)=x;
    dP(2,:)=ones(1,n);
    for i=2:order
        P(i+1,:)=((2*i-1)/i)*x.*P(i,:)-((i-1)/i)*P(i-1,:);
        dP(i+1,:)=((2*i-1)/i)*P(i,:)+((2*i-1)/i)*x.*dP(i,:)-((i-1)/i)*dP(i-1,:);
    end
end
%val=P(order+1,:);
%dval=dP(order+1,:);
val=P;
dval=dP;
end