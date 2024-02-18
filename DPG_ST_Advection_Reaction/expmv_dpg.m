function [f,s,m,mv,mvd,unA] = ...
          expmv_dpg(t,A1,A2,A3,A4,W,J,b,M,prec,full_term,prnt)
%EXPMV   Matrix exponential times vector or matrix.
%   [F,S,M,MV,MVD] = EXPMV(t,A,B,[],PREC) computes EXPM(t*A)*B without
%   explicitly forming EXPM(t*A). PREC is the required accuracy, 'double',
%   'single' or 'half', and defaults to CLASS(A).
%   A total of MV products with A or A^* are used, of which MVD are
%   for norm estimation.
%   The full syntax is
%     [f,s,m,mv,mvd,unA] = expmv(t,A,b,M,prec,shift,bal,full_term,prnt).
%   unA = 1 if the alpha_p were used instead of norm(A).
%   If repeated invocation of EXPMV is required for several values of t
%   or B, it is recommended to provide M as an external parameter as
%   M = SELECT_TAYLOR_DEGREE(A,b,m_max,p_max,prec,shift,bal,true).
%   This also allows choosing different m_max and p_max.

%   Reference: A. H. Al-Mohy and N. J. Higham. Computing the action of the
%   matrix exponential, with an application to exponential integrators.
%   SIAM J. Sci. Comput., 33(2):488--511, 2011.  Algorithm 3.2.

%   Awad H. Al-Mohy and Nicholas J. Higham, November 9, 2010.

if nargin < 12 || isempty(prnt), prnt = false; end
if nargin < 11 || isempty(full_term), full_term = false; end
if nargin < 10 || isempty(prec), prec = class(A1); end
if nargin < 9 || isempty(M)
   tt = 1;
   [M,mvd,alpha,unA] = select_taylor_degree_dpg(t*A1,t*A2,A3,A4,t*W,t*J,b,[],[],prec);
   mv = mvd;
else
   tt = t; mv = 0; mvd = 0;
end

switch prec
    case 'double', tol = 2^(-53);
    case 'single', tol = 2^(-24);
    case 'half',   tol = 2^(-10);
end

n=length(A1);

s = 1;
if t == 0
    m = 0;
else
    [m_max,p] = size(M);
     U = diag(1:m_max);
     C = ( (ceil(abs(tt)*M))'*U );
     C (C == 0) = inf;
     if p > 1
         [cost m] = min(min(C)); % cost is the overall cost.
     else
         [cost m] = min(C);  % when C is one column. Happens if p_max = 2.
     end
     if cost == inf; cost = 0; end
     s = max(cost/m,1);
end
f = b;
if prnt, fprintf('m = %2.0f, s = %g, m_actual = ', m, s), end
for i = 1:s
    c1 = norm(b,inf);
    for k = 1:m

        b = (t/(s*k))*actionAb(A1,A2,A3,A4,W,J,b,n);
        mv = mv + 4;
        f =  f + b;
        c2 = norm(b,inf);
        if ~full_term
%            if prnt, fprintf(' %9.2e, \n', (c1+c2)/norm(f,inf)), pause, end
            if c1 + c2 <= tol*norm(f,inf)
                if prnt, fprintf(' %2.0f, ', k), end
                break
            end
            c1 = c2;
        end

    end
    b = f;
end
if prnt, fprintf('\n'), end
end

function B=actionAb(A1,A2,A3,A4,W,J,b,n)
    B0 = A4*b(1:n,:); 
    B1 = A3*B0;
    B=[A1*b(1:n,:)-A2*B1+W*b(n+1:end,:);J*b(n+1:end,:)];
end
