%% Programa para obtener K y f correspondientes al elemento de barra de 
%% 2 nodos unidimensional, utilizando funciones de forma locales

syms x x1 x2 E A L b;          % definicion de las variables simbolicas
 
% Defino las funciones de forma
x2 = x1 + L;
N1 = (x2-x)/L;                 N2 = (x-x1)/L;

N = [N1 N2];                   % matriz de funciones de forma
B = [diff(N1,x) diff(N2,x)];   % matriz de deformacion
D = E*A;                       % matriz constitutiva

% Matriz de rigidez (ecuacion 2.76b)
K = int(B.'*D*B, x, x1, x2);
disp('K = '); pretty(K)
 
% Vector de fuerzas nodales equivalentes (ecuacion 2.76b)
f = int(N.'*b, x, x1, x2);
disp('f = '); pretty(f)
