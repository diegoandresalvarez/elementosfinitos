%% Programa para obtener K y f correspondientes al elemento de barra de 
%% 2 nodos unidimensional, utilizando funciones de forma globales

syms x L E A b N1(x) N2(x) N3(x) N4(x)

% Defino las funciones de forma
N = [N1(x), N2(x), N3(x), N4(x)]; % matriz de funciones de forma
B = diff(N,x);                    % matriz de deformacion
D = E*A;                          % matriz constitutiva

% Matriz de rigidez (ecuacion 2.83)
K = simplify(int(B.'*D*B, x, 0, L))
disp('K = '); pretty(K);

% Vector de fuerzas nodales equivalentes (ecuacion 2.83)
f = simplify(int(N.'*b, x, 0, L))
disp('f = '); pretty(f);
