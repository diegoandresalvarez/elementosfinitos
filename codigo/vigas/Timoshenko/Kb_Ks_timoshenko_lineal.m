% Calculo de las funciones de forma del elemento de viga de Timoshenko
% lineal

%% Borro memoria, pantalla y cierro gráficos
clear, clc, close all

%% Definición de variables
syms x xi E I G Aast L w1 t1 w2 t2 fz m

%% Vector de movimientos nodales
a = [w1; t1; w2; t2];

%% Defino las funciones de forma
N1 = (1-xi)/2;
N2 = (1+xi)/2;

%% Se define el inverso del Jacobiano y el Jacobiano
dx_dxi = L/2;
dxi_dx = 2/L;

%% Defino las matrices de deformacion
Bb = [                  0, dxi_dx*diff(N1,xi),                  0, dxi_dx*diff(N2,xi) ]
Bs = [ dxi_dx*diff(N1,xi),                -N1, dxi_dx*diff(N2,xi),                -N2 ]

%% Defino las matrices de rigidez
Kb = int(Bb.'*E*I*Bb*dx_dxi,   xi,-1,1);
Ks = int(Bs.'*G*Aast*Bs*dx_dxi,xi,-1,1);
Us_ex = simplify(a.'*Ks*a/2);

disp('Solucion exacta = ')
disp('Kb = (E*I/L) * '),    pretty(Kb/(E*I/L))
disp('Ks = (G*Aast/L) * '), pretty(Ks/(G*Aast/L))

%% Evaluo la curvatura
kappa = simplify(dxi_dx*(diff(N1,xi)*t1 + diff(N2,xi)*t2));
disp('kappa = '), pretty(kappa)

%% Evaluo gamma_xz
gxz =  expand(dxi_dx*(diff(N1,xi)*w1 + diff(N2,xi)*w2) - (N1*t1 + N2*t2));
% Us_ex = int(G*Aast*gxz^2*L/2,xi,-1,1)/2
disp('gxz = '), pretty(collect(gxz,xi))
fprintf('\n\n\n\n\n\n');

%% Evaluo las matrices con una cuadratura de Gauss-Legendre de orden 1
w1 = 2; xi1 = 0;
Kb = subs(Bb.'*E*I*Bb*dx_dxi, xi, xi1)*w1;
Ks = subs(Bs.'*G*Aast*Bs*dx_dxi, xi, xi1)*w1;
Us_1 = simplify(a.'*Ks*a/2);

disp('Integrando con GL de orden 1 = ')
disp('Kb = (E*I/L) * '),    pretty(Kb/(E*I/L))
disp('Ks = (G*Aast/L) * '), pretty(Ks/(G*Aast/L))
fprintf('\n\n\n\n\n\n');

%% Evaluo las matrices con una cuadratura de Gauss-Legendre de orden 2
w1 = 1; xi1 = -sym(sqrt(1/3));
w2 = 1; xi2 = +sym(sqrt(1/3));
Kb = simplify(subs(Bb.'*E*I*Bb*dx_dxi,   xi, xi1)*w1 + ...
              subs(Bb.'*E*I*Bb*dx_dxi,   xi, xi2)*w2);
           
Ks = simplify(subs(Bs.'*G*Aast*Bs*dx_dxi,xi, xi1)*w1 + ...
              subs(Bs.'*G*Aast*Bs*dx_dxi,xi, xi2)*w2);
Us_2 = simplify(a.'*Ks*a/2);          

disp('Integrando con GL de orden 2 = ')
disp('Kb = (E*I/L) * '),    pretty(Kb/(E*I/L))
disp('Ks = (G*Aast/L) * '), pretty(Ks/(G*Aast/L))

%% Calculo el vector de fuerzas nodales equivalentes correspondientes a una
%% carga distribuida de magnitud q constante
f1 = int(N1*[fz; m]*dx_dxi, xi, -1, 1);
f2 = int(N2*[fz; m]*dx_dxi, xi, -1, 1);
fe = [f1; f2];
disp('fe = '), pretty(fe)

%% Calculo el vector de fuerzas nodales equivalentes correspondientes a una
%% carga distribuida trapezoidal y sin momentos distribuidos
syms q1 q2
fz = N1*q1 + N2*q2;
m = 0;
f1 = int(N1*[fz; m]*dx_dxi, xi, -1, 1);
f2 = int(N2*[fz; m]*dx_dxi, xi, -1, 1);
fe = [f1; f2];
disp('fe = '), pretty(fe)

%% Bye, bye!!
return