%%%Convergence rate%%%
% T0=readtable('ConvergenceFieldsTime2.txt');
% Tmat0=table2array(T0);
% x0=log(Tmat0(2:end,2)./Tmat0(1:end-1,2));
% y0=log(Tmat0(2:end,1)./Tmat0(1:end-1,1));
% ConvRates0=x0./y0
% 
T0=readtable('ConvergenceTraceL2Space2.txt');
Tmat0=table2array(T0);
x0=log(Tmat0(2:end,2)./Tmat0(1:end-1,2));
y0=log(Tmat0(2:end,1)./Tmat0(1:end-1,1));
ConvRates0=x0./y0

% T0=readtable('ConvergenceTraceMaxSpace2.txt');
% Tmat0=table2array(T0);
% x0=log(Tmat0(2:end,2)./Tmat0(1:end-1,2));
% y0=log(Tmat0(2:end,1)./Tmat0(1:end-1,1));
% ConvRates0=x0./y0
