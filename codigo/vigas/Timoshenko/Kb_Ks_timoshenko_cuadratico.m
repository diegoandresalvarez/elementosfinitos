clear, clc, close all

%% Definicion de variables
syms xi x1 x2 x3 L E I Aast G
syms w1 w2 w3 t1 t2 t3

%% Vector de movimientos nodales
a = [w1; t1; w2; t2; w3; t3];

%% Defino las posiciones de los nodos
x3 = x1+L;
x2 = (x1+x3)/2;

%% Funciones de forma Lagrangianas
N1 = calc_N([-1 0 1], [1 0 0], xi); % = xi*(xi-1)/2
N2 = calc_N([-1 0 1], [0 1 0], xi); % = (1+xi)*(1-xi)
N3 = calc_N([-1 0 1], [0 0 1], xi); % = xi*(xi+1)/2

%% Interpolacion de la geometria, la flecha y el giro
x = simplify(N1*x1 + N2*x2 + N3*x3);
%w = simplify(N1*w1 + N2*w2 + N3*w3);
%t = simplify(N1*t1 + N2*t2 + N3*t3);

%% Derivadas
dx_dxi = diff(x,xi);
syms xx % creo una variable "dummy" para despejar xi de x = N1*x1 + N2*x2 + N3*x3;
dxi_dx = diff(solve(xx == x,xi),xx);

%% Defino las matrices de deformacion
Bb = [ 0, diff(N1,xi)*dxi_dx  , 0, diff(N2,xi)*dxi_dx  , 0, diff(N3,xi)*dxi_dx   ];
Bs = [ diff(N1,xi)*dxi_dx, -N1, diff(N2,xi)*dxi_dx, -N2, diff(N3,xi)*dxi_dx, -N3 ];

disp('Bb = (2/L) * '),    pretty(simplify(Bb/(2/L)))
disp('Bs = (2/L) * '),    pretty(simplify(Bs/(2/L)))

%% Integro exactamente las matrices de rigidez
disp('Integral exacta de las matrices de rigidez = ');
Kb = int(Bb.'*E*I*Bb*L/2,   xi,-1,1);
Ks = int(Bs.'*G*Aast*Bs*L/2,xi,-1,1);
Us_ex = simplify(a.'*Ks*a/2);

disp('Kb = ((E*I)/(3*L)) * '),    pretty(Kb/(E*I/(3*L)))
disp('Ks = ((G*Aast)/(9*L)) * '), pretty(Ks/(G*Aast/(9*L)))

%% Evaluo la curvatura
kappa = simplify(dxi_dx*(diff(N1,xi)*t1 + diff(N2,xi)*t2 + + diff(N3,xi)*t3));
disp('kappa = '), pretty(kappa)

%% Evaluo gamma_xz
gxz =  simplify(dxi_dx*(diff(N1,xi)*w1 + diff(N2,xi)*w2 + diff(N3,xi)*w3) - (N1*t1 + N2*t2 + + N3*t3));
% Us_ex = int(G*Aast*gxz^2*L/2,xi,-1,1)/2
disp('gxz = '), pretty(collect(gxz,xi))
fprintf('\n\n\n\n\n\n');

%% Integro las matrices con una cuadratura de Gauss-Legendre de orden 2
disp('Integral con una cuadratura de Gauss-Legendre de orden 2 = ');
xi1 = -sym(sqrt(1/3));   w1 = 1;
xi2 = +sym(sqrt(1/3));   w2 = 1;
Kb = simplify(subs(Bb.'*E*I*Bb*L/2,   xi,xi1)*w1 + ...
              subs(Bb.'*E*I*Bb*L/2,   xi,xi2)*w2);
         
Ks = simplify(subs(Bs.'*G*Aast*Bs*L/2,xi,xi1)*w1 + ...
              subs(Bs.'*G*Aast*Bs*L/2,xi,xi2)*w2);
Us_2 = simplify(a.'*Ks*a/2);          

disp('Kb = ((E*I)/(3*L)) * '),    pretty(Kb/(E*I/(3*L)))
disp('Ks = ((G*Aast)/(9*L)) * '), pretty(Ks/(G*Aast/(9*L)))

%% Integro las matrices con una cuadratura de Gauss-Legendre de orden 3
disp('Integral con una cuadratura de Gauss-Legendre de orden 3 = ');
xi1 = -sym(sqrt(3/5));   w1 = sym(5/9);
xi2 =  sym(0);           w2 = sym(8/9);
xi3 = +sym(sqrt(3/5));   w3 = sym(5/9);
Kb = simplify(subs(Bb.'*E*I*Bb*L/2,   xi,xi1)*w1 + ...
              subs(Bb.'*E*I*Bb*L/2,   xi,xi2)*w2 + ...
              subs(Bb.'*E*I*Bb*L/2,   xi,xi3)*w3) ;
         
Ks = simplify(subs(Bs.'*G*Aast*Bs*L/2,xi,xi1)*w1 + ...
              subs(Bs.'*G*Aast*Bs*L/2,xi,xi2)*w2 + ...            
              subs(Bs.'*G*Aast*Bs*L/2,xi,xi3)*w3);
Us_3 = simplify(a.'*Ks*a/2);          

disp('Kb = ((E*I)/(3*L)) * '),    pretty(Kb/(E*I/(3*L)))
disp('Ks = ((G*Aast)/(9*L)) * '), pretty(Ks/(G*Aast/(9*L)))

%% Bye, bye!!!
return

%% -------------------------------------------------------------------------
%% Calcular correctamente los polinomios de las funciones de forma 1D
function N = calc_N(xp, yp, var)
    % se ve verifican los tamanios de los vectores xp y yp
    nx = length(xp);
    ny = length(yp);
    assert(nx == ny, 'Los vectores xp y yp deben tener el mismo tamanio');

    % se calculan los coeficientes de los polinomios
    c = polyfit(xp, yp, nx-1);
    
    % se eliminan los errores en la aproximacion numerica, haciendo los
    % coeficientes demasiado pequenios igual a cero
    c(abs(c) < 1e-10) = 0;
    
    % con los coeficientes corregidos se calculan las funciones de forma
    N = poly2sym(c, var);
end
