function  [M,mv,alpha,unA] = ...
           select_taylor_degree_dpg(A1,A2,A3,A4,W,J,b,m_max,p_max,prec,force_estm)
%SELECT_TAYLOR_DEGREE   Select degree of Taylor approximation.
%   [M,MV,alpha,unA] = SELECT_TAYLOR_DEGREE(A,b,m_max,p_max) forms a matrix M
%   for use in determining the truncated Taylor series degree in EXPMV based 
%   on parameters m_max and p_max.
%   MV is the number of matrix-vector products with A or A^* computed.

if nargin < 11, force_estm = false; end
if nargin < 9 || isempty(p_max), p_max = 8; end
if nargin < 8 || isempty(m_max), m_max = 55; end

if p_max < 2 || m_max > 60 || m_max + 1 < p_max*(p_max - 1)
    error('>>> Invalid p_max or m_max.')
end

if nargin < 10 || isempty(prec), prec = class(A1); end
switch prec
    case 'double'
        load theta_taylor
    case 'single'
        load theta_taylor_single
    case 'half'
        load theta_taylor_half
end

mv = 0;
if ~force_estm, normA = normAm_dpg(A1,A2,A3,A4,W,J,1); end

if ~force_estm && normA <= 4*theta(m_max)*p_max*(p_max + 3)/(m_max*size(b,2))
    % Base choice of m on normA, not the alpha_p.
    unA = 1;
    c = normA;
    alpha = c*ones(p_max-1,1);
 else
    unA = 0;
    eta = zeros(p_max,1); alpha = zeros(p_max-1,1);
    for p = 1:p_max
        [c,k] = normAm_dpg(A1,A2,A3,A4,W,J,p+1);
        c = c^(1/(p+1));
        mv = mv + k;
        eta(p) = c;
    end
    for p = 1:p_max-1
        alpha(p) = max(eta(p),eta(p+1));
    end
end
M = zeros(m_max,p_max-1);
for p = 2:p_max
    for m = p*(p-1)-1 : m_max
        M(m,p-1) = alpha(p-1)/theta(m);
    end
end
