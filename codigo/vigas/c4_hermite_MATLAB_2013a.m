%% Interpolacion polinomica de Hermite
clear, clc

%% Definicion de variables
syms xi u1 u2 u1p u2p
x1 = -1; x2 =  1;

%% Funciones de forma Lagrangianas
L1 = poly2sym(polyfit([-1 1],[1 0],1),xi);
L2 = poly2sym(polyfit([-1 1],[0 1],1),xi);

%% Se calculan los polinomios Pi y Qi
P1 = (1 - 2*subs(diff(L1,xi),xi,x1)*(xi-x1))*L1^2;
P2 = (1 - 2*subs(diff(L2,xi),xi,x2)*(xi-x2))*L2^2;
Q1 = (xi-x1)*L1^2; 
Q2 = (xi-x2)*L2^2;

%% Se calcula el polinomio interpolador
u = simple(P1*u1 + P2*u2 + Q1*u1p + Q2*u2p);

%% Se calculan las funciones de forma
N1  = expand(subs(u, {u1, u1p, u2, u2p} , {1, 0, 0, 0}));
N1b = expand(subs(u, {u1, u1p, u2, u2p} , {0, 1, 0, 0}));
N2  = expand(subs(u, {u1, u1p, u2, u2p} , {0, 0, 1, 0}));
N2b = expand(subs(u, {u1, u1p, u2, u2p} , {0, 0, 0, 1}));

%% Se muestran los resultados
disp('N_1(xi)  = '); pretty(N1);    
disp('Nb_1(xi) = '); pretty(N1b);
disp('N_2(xi)  = '); pretty(N2);    
disp('Nb_2(xi) = '); pretty(N2b);

%% bye bye