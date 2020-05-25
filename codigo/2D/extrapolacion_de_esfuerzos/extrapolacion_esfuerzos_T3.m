clear, clc

%% Se define el elemento finito
% Numeracion local del EF triangular de 3 nodos:
%  ^ eta
%  |
%  3
%  |\
%  | \ 
%  |  \
%  |   \
%  1----2---> xi

% Coordenadas de los nodos
xnod = [ ...
%  xi   eta     % nodo   
    0    0      %  1
    1    0      %  2
    0    1  ];  %  3

%% se define el orden de la cuadratura
orden = 3;
GP = TriGaussPoints(orden); nGP = size(GP,1);
xi_gl = GP(:,1);   eta_gl = GP(:,2);

%% se calcula la matriz de interpolacion A
A1 = sym([ ones(nGP, 1) xi_gl eta_gl xi_gl.*eta_gl ]);
A2 = sym([ ones(3,1) xnod(:,1)  xnod(:,2)  xnod(:,1).*xnod(:,2) ]);
A = simplify(A2/A1); %= A2*inv(A1)

%% y se imprime
disp('La matriz de interpolacion es:')
disp(A)