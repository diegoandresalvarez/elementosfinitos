% Calculo de las funciones de forma del elemento de viga de Timoshenko
% lineal

%%
clear, clc, close all

syms x xi E I G Aast L w1 t1 w2 t2 q

%% Defino las funciones de forma
N1 = (1-xi)/2;
N2 = (1+xi)/2;

%% Defino las matrices de deformacion
Bf = [ 0,                 (2/L)*diff(N1,xi), 0,                 (2/L)*diff(N2,xi) ]
Bc = [ (2/L)*diff(N1,xi), -N1,               (2/L)*diff(N2,xi), -N2               ]

%% Defino las matrices de rigidez
Kf = int(Bf.'*E*I*Bf*L/2,   xi,-1,1);
Kc = int(Bc.'*G*Aast*Bc*L/2,xi,-1,1);

disp('Solucion exacta = ')
disp('Kf = (E*I/L) * '),    pretty(Kf/(E*I/L))
disp('Kc = (G*Aast/L) * '), pretty(Kc/(G*Aast/L))
fprintf('\n\n\n\n\n\n');

%% Evaluo las matrices con una cuadratura de Gauss-Legendre de orden 1
w = 2; % para xi = 0
% El comando subs no funciona, de todos modos no lo necesitamos:
Kf = (Bf.'*E*I*Bf*L/2   )*w; % Kf = subs(Bf.'*E*I*Bf*L/2,   xi, 0)*w;
Kc = subs(Bc.'*G*Aast*Bc*L/2,xi, 0)*w;

disp('Integrando con GL de orden 1 = ')
disp('Kf = (E*I/L) * '),    pretty(Kf/(E*I/L))
disp('Kc = (G*Aast/L) * '), pretty(Kc/(G*Aast/L))
fprintf('\n\n\n\n\n\n');

%% Evaluo las matrices con una cuadratura de Gauss-Legendre de orden 2
w = 1; % para xi = -sqrt(1/3) y sqrt(1/3)
Kf = simple(subs(Bf.'*E*I*Bf*L/2,   xi,-sym(sqrt(1/3)))*w + subs(Bf.'*E*I*Bf*L/2,   xi,sym(sqrt(1/3)))*w);
Kc = simple(subs(Bc.'*G*Aast*Bc*L/2,xi,-sym(sqrt(1/3)))*w + subs(Bc.'*G*Aast*Bc*L/2,xi,sym(sqrt(1/3)))*w);

disp('Integrando con GL de orden 2 = ')
disp('Kf = (E*I/L) * '),    pretty(Kf/(E*I/L))
disp('Kc = (G*Aast/L) * '), pretty(Kc/(G*Aast/L))

%% Evaluo la curvatura
chi = simple((2/L)*(diff(N1,xi)*t1 + diff(N2,xi)*t2));
disp('CHI = '), pretty(chi)

%% Evaluo gamma_xz
gxz =  simple((2/L)*(diff(N1,xi)*w1 + diff(N2,xi)*w2) - (N1*t1 + N2*t2));
disp('gxz = '), pretty(gxz)

%% Calculo el vector de fuerzas nodales equivalentes correspondientes a una
%% carga distribuida de magnitud q constante
Nw = [ N1  0  N2  0 ];

fe = int(Nw.'*q*L/2, xi,-1,1)

return %bye, bye!
