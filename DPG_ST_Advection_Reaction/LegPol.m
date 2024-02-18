function val = LegPol(x,order)
P=zeros(order,size(x,2));
if order==1
    P(1,:)=ones(1,size(x,2));
elseif order==2
    P(2,:)=x;
else 
    P(1,:)=ones(1,size(x,2));
    P(2,:)=x;
    for i=3:order
        n=i-1;
        P(i,:)=((2*n-1)/n)*x.*P(i-1,:)-((n-1)/n)*P(i-2,:);
    end
end
val=P(order,:);
end