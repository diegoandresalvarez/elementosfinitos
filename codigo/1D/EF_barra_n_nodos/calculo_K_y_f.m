%% Programa para obtener K y f correspondientes al elemento isoparametrico 
%% lagrangiano cuadratico de tres nodos unidimensional 

%% Definicion de variables
syms xi x1 x2 x3 L E A b

%% Defino las posiciones de los nodos
x3 = x1+L;
x2 = (x1+x3)/2;

%% Funciones de forma lagrangianas
%N1 = poly2sym(polyfit([-1 0 1],[1 0 0],2),xi);   % = xi*(xi-1)/2;
%N2 = poly2sym(polyfit([-1 0 1],[0 1 0],2),xi);   % = (1+xi)*(1-xi);
%N3 = poly2sym(polyfit([-1 0 1],[0 0 1],2),xi);   % = xi*(xi+1)/2;

% Se calculan los coeficientes de los polinomios
c1 = polyfit([-1 0 1],[1 0 0],2);
c2 = polyfit([-1 0 1],[0 1 0],2);
c3 = polyfit([-1 0 1],[0 0 1],2);

% Se eliminan los errores en la aproximacion numerica, haciendo los
% coeficientes demasiado pequenios igual a cero
c1(abs(c1) < 1e-10) = 0;
c2(abs(c2) < 1e-10) = 0;
c3(abs(c3) < 1e-10) = 0;

% con los coeficientes corregidos se calculan las funciones de forma
N1 = poly2sym(c1, xi);
N2 = poly2sym(c2, xi);
N3 = poly2sym(c3, xi);

%% Interpolacion de la geometria y sus derivadas
x   = simplify(N1*x1 + N2*x2 + N3*x3);       
dx_dxi = diff(x, xi);

% NOTA: el siguiente comando solo se puede realizar si x(xi) es una función
% continua e inyectiva. Este es nuestro caso:
% https://en.wikipedia.org/wiki/Inverse_functions_and_differentiation
dxi_dx = 1/dx_dxi;
% recuerde que se debe garantizar que dx_dxi>0 y dxi_dx>0

%% Definicion de la matriz de forma N y matriz de deformacion del elemento B
N = [N1 N2 N3];
% B = simple([diff(N1,xi) diff(N2,xi) diff(N3,xi)]*dxi_dx);
B = diff(N,xi)*dxi_dx;

%% "matriz constitutiva"
D = E*A;

%% Calculo la matriz de rigidez del elemento
K = int(B.'*D*B*dx_dxi, xi, -1, +1);
disp('K = ((A*E)/(3*L))*')
pretty(K/((A*E)/(3*L)))

%% Calculo la matriz de fuerzas nodales equivalentes del elemento
f = int(N.'*b*dx_dxi, xi, -1, +1);
disp('f = ((b*L)/6)*')
pretty(f/((b*L)/6))
