function [c,mv] = normAm_dpg(A1,A2,A3,A4,W,J,m)
%NORMAM   Estimate of 1-norm of power of matrix.
%   NORMAM_DPG(A,m) estimates norm(A^m,1) where A=A1-A2*A3*A4

%   [C,MV] = NORMAM_DPG(A,m) returns the estimate C and the number MV of
%   matrix-vector products computed involving A or A^*.

t = 1; % Number of columns used by NORMEST1.

n = length(A1);
q = length(J);

[c,v,w,it] = normest1(@afun_power,t);
mv = it(2)*t*m*4;

  function Z = afun_power(flag,X)
       %AFUN_POWER  Function to evaluate matrix products needed by NORMEST1.

       if isequal(flag,'dim')
          Z = n+q;
       elseif isequal(flag,'real')
          Z = isreal(A1);
       else
          p = size(X,1);
          if p ~= n+q, error('Dimension mismatch'), end

          if isequal(flag,'notransp')
             for i = 1:m 
                 X0 = A4*X(1:n,:); 
                 X1 = A3*X0;
                 X=[A1*X(1:n,:)-A2*X1+W*X(n+1:end,:);J*X(n+1:end,:)];
             end
          elseif isequal(flag,'transp')
             for i = 1:m
                 X0 = A2'*X(1:n,:); 
                 X1 = A3'*X0;
                 X=[A1'*X(1:n,:)-A4'*X1;W'*X(1:n,:)+J'*X(n+1:end,:)];
             end
          end

          Z = X;

       end

  end
end
