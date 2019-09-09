%% Programa para obtener K y f correspondientes al elemento isoparametrico 
%% lagrangiano cuadratico (de tres nodos) unidimensional 

clear all; clc

%% Definicion de variables
syms xi x1 x2 x3 L E A b

%% Defino las posiciones de los nodos
x3 = x1+L;
x2 = (x1+x3)/2;

%% Funciones de forma Lagrangianas
N1 = poly2sym(polyfit([-1 0 1],[1 0 0],2),xi);   % = xi*(xi-1)/2;
N2 = poly2sym(polyfit([-1 0 1],[0 1 0],2),xi);   % = (1+xi)*(1-xi);
N3 = poly2sym(polyfit([-1 0 1],[0 0 1],2),xi);   % = xi*(xi+1)/2;

%% Interpolacion de la geometria y sus derivadas
x   = simple(N1*x1 + N2*x2 + N3*x3);       dx_dxi = diff(x,   xi);
xi_ = solve(['x = ' char(x)],xi);          dxi_dx = diff(xi_, 'x');
% recuerde que se debe garantizar que dx_dxi>0 y dxi_dx>0

%% Definicion de la matriz de forma N y matriz de deformacion del elemento B
N = [N1 N2 N3];
B = simple([diff(N1,xi) diff(N2,xi) diff(N3,xi)]*dxi_dx);

%% "matriz constitutiva"
D = E*A;

%% Calculo la matriz de rigidez del elemento
K = int(B.'*D*B*dx_dxi,xi,-1,1);
disp('K = ((A*E)/(6*L))*')
pretty(K/((A*E)/(6*L)))

%% Calculo la matriz de fuerzas nodales equivalentes del elemento
f = int(N.'*b*dx_dxi,xi,-1,1);
disp('f = ((b*L)/6)*')
pretty(f/((b*L)/6))

%% Tarea
% hacer programa para obtener K y f correspondientes al elemento
% isoparametrico lagrangiano cúbico (de cuatro nodos) unidimensional 
