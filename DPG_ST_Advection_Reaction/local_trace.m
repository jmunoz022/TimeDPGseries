function l=local_trace(order)
%%%Function that stores the local gramm matrix for each element%%% 

%xcoord=x coordinates of the element
%nquad=number of quadrature points

%%%Integration%%%
l=zeros(order+1,2);
l(:,1)=-Shape1Q(-1,order);
l(:,2)=Shape1Q(1,order);

end

