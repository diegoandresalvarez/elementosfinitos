clear, clc, close all

%% Definicion de variables
syms xi x1 x2 x3 L E I Aast G
% syms w1 w2 w3 t1 t2 t3

%% Defino las posiciones de los nodos
x3 = x1+L;
x2 = (x1+x3)/2;

%% Funciones de forma Lagrangianas
N1 = poly2sym(polyfit([-1 0 1],[1 0 0],2),xi);   % = xi*(xi-1)/2;
N2 = poly2sym(polyfit([-1 0 1],[0 1 0],2),xi);   % = (1+xi)*(1-xi);
N3 = poly2sym(polyfit([-1 0 1],[0 0 1],2),xi);   % = xi*(xi+1)/2;

%% Interpolacion de la geometria, la flecha y el giro
x = simplify(N1*x1 + N2*x2 + N3*x3);
%w = simplify(N1*w1 + N2*w2 + N3*w3);
%t = simplify(N1*t1 + N2*t2 + N3*t3);

%% Derivadas
dx_dxi = diff(x,xi);
syms xx % creo una variable "dummy" para despejar xi de x = N1*x1 + N2*x2 + N3*x3;
dxi_dx = diff(solve(xx - x,xi),xx);

%% Defino las matrices de deformacion
Bf = [ 0, diff(N1,xi)*dxi_dx  , 0, diff(N2,xi)*dxi_dx  , 0, diff(N3,xi)*dxi_dx   ];
Bc = [ diff(N1,xi)*dxi_dx, -N1, diff(N2,xi)*dxi_dx, -N2, diff(N3,xi)*dxi_dx, -N3 ];

disp('Bf = (2/L) * '),    pretty(simplify(Bf/(2/L)))
disp('Bc = (2/L) * '),    pretty(simplify(Bc/(2/L)))

%% Integro exactamente las matrices de rigidez
disp('Integral exacta de las matrices de rigidez = ');
Kf = int(Bf.'*E*I*Bf*L/2,   xi,-1,1);
Kc = int(Bc.'*G*Aast*Bc*L/2,xi,-1,1);

disp('Kf = ((E*I)/(3*L)) * '),    pretty(Kf/(E*I/(3*L)))
disp('Kc = ((G*Aast)/(9*L)) * '), pretty(Kc/(G*Aast/(9*L)))

%% Integro las matrices con una cuadratura de Gauss-Legendre de orden 2
disp('Integral con una cuadratura de Gauss-Legendre de orden 2 = ');
xi1 = -sym(sqrt(1/3));   w1 = 1;
xi2 = +sym(sqrt(1/3));   w2 = 1;
Kf = simplify(subs(Bf.'*E*I*Bf*L/2,   xi,xi1)*w1 + ...
              subs(Bf.'*E*I*Bf*L/2,   xi,xi2)*w2) ;
         
Kc = simplify(subs(Bc.'*G*Aast*Bc*L/2,xi,xi1)*w1 + ...
              subs(Bc.'*G*Aast*Bc*L/2,xi,xi2)*w2);

disp('Kf = ((E*I)/(3*L)) * '),    pretty(Kf/(E*I/(3*L)))
disp('Kc = ((G*Aast)/(9*L)) * '), pretty(Kc/(G*Aast/(9*L)))

%% Integro las matrices con una cuadratura de Gauss-Legendre de orden 3
disp('Integral con una cuadratura de Gauss-Legendre de orden 3 = w1');
xi1 = -sym(sqrt(3/5));   w1 = sym(5/9);
xi2 =  sym(0);           w2 = sym(8/9);
xi3 = +sym(sqrt(3/5));   w3 = sym(5/9);
Kf = simplify(subs(Bf.'*E*I*Bf*L/2,   xi,xi1)*w1 + ...
              subs(Bf.'*E*I*Bf*L/2,   xi,xi2)*w2 + ...
              subs(Bf.'*E*I*Bf*L/2,   xi,xi3)*w3) ;
         
Kc = simplify(subs(Bc.'*G*Aast*Bc*L/2,xi,xi1)*w1 + ...
              subs(Bc.'*G*Aast*Bc*L/2,xi,xi2)*w2 + ...            
              subs(Bc.'*G*Aast*Bc*L/2,xi,xi3)*w3);

disp('Kf = ((E*I)/(3*L)) * '),    pretty(Kf/(E*I/(3*L)))
disp('Kc = ((G*Aast)/(9*L)) * '), pretty(Kc/(G*Aast/(9*L)))

return %bye, bye!
