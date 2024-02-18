x=linspace(-1,1,100);

[P0,dP0]=Shape1Q(x,0);
[P1,dP1]=Shape1Q(x,1);
[P2,dP2]=Shape1Q(x,2);
[P3,dP3]=Shape1Q(x,3);
[P4,dP4]=Shape1Q(x,4);


[L1,dL1]=Shape1H(x,1);
[L2,dL2]=Shape1H(x,2);
[L3,dL3]=Shape1H(x,3);
[L4,dL4]=Shape1H(x,4);

figure
plot(x,P0,x,P1,x,P2,x,P3,x,P4)
title('Legendre')

figure
plot(x,dP0,x,dP1,x,dP2,x,dP3,x,dP4)
title('Legendre derivative')


figure
plot(x,L1(1,:),x,L1(2,:),x,L2,x,L3,x,L4)
title('Lobatto')

figure
plot(x,dL1(1,:),x,dL1(2,:),x,dL2,x,dL3,x,dL4)
title('Lobatto derivative')



