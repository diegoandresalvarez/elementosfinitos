% Ejemplo 4.13 Oñate (1995)
% Ejemplo 2.9 Oñate (2013)

%% A partir de la interpolacion
% w = pol grado 2
% t = pol grado 2
% e imponiendo que gxz varie linealmente entre los valores
% gxz_1 = gxz(xi=-1/sqrt(3)) y 
% gxz_2 = gxz(xi=+1/sqrt(3))
% obtenga un elemento de viga de Timoshenko tres nodos

clear, clc

syms xi w1 w2 w3 t1 t2 t3 E I Aast G L
dx_dxi = L/2;
dxi_dx = 2/L;

% -1                              0                              1 
%  x------------------------------x------------------------------x--> xi
%  w1                             w2                             w3
%  t1                             t2                             t3

%% Funciones de forma Lagrangianas
N1_2 = calc_N([-1 0 1], [1 0 0], xi); % = xi*(xi-1)/2
N2_2 = calc_N([-1 0 1], [0 1 0], xi); % = (1+xi)*(1-xi)
N3_2 = calc_N([-1 0 1], [0 0 1], xi); % = xi*(xi+1)/2

%% Se definen w y t
w = N1_2*w1 + N2_2*w2 + N3_2*w3;
t = N1_2*t1 + N2_2*t2 + N3_2*t3;

%% Se define el vector de movimientos nodales del elemento ae
ae = {w1,t1,w2,t2,w3,t3};

%% Se recalcula gxz y se calcula la matriz Bs
gxz = expand(diff(w,xi)*dxi_dx - t);
Bs = simplify([ subs(gxz,ae,{1,0,0,0,0,0}), ...
                subs(gxz,ae,{0,1,0,0,0,0}), ...
                subs(gxz,ae,{0,0,1,0,0,0}), ...
                subs(gxz,ae,{0,0,0,1,0,0}), ...
                subs(gxz,ae,{0,0,0,0,1,0}), ...
                subs(gxz,ae,{0,0,0,0,0,1}) ]);

%a = sym(1/sqrt(3))
syms a
%gxz_1 = subs(Bs,xi,-a)*[w1;t1;w2;t2;w3;t3]
%gxz_2 = subs(Bs,xi,+a)*[w1;t1;w2;t2;w3;t3]

% Los puntos de evaluacion seran xi1 = -1/sqrt(3) y xi2 = -1/sqrt(3)
gxz_1 = subs(gxz,xi,-a)
gxz_2 = subs(gxz,xi,+a)

Ng1 = poly2sym(polyfit([-1/sqrt(3) 1/sqrt(3)],[1 0],1),xi); Ng1 = subs(Ng1,3^(1/2),1/a);
Ng2 = poly2sym(polyfit([-1/sqrt(3) 1/sqrt(3)],[0 1],1),xi); Ng2 = subs(Ng2,3^(1/2),1/a);
% toca hacer ese truco, porque o sino el toolbox simbolico del MATLAB no
% produce buenos resultados

gxz_lin = Ng1*gxz_1 + Ng2*gxz_2;
Bs_sustitutiva = simplify([ subs(gxz_lin,ae,{1,0,0,0,0,0}), ...
                            subs(gxz_lin,ae,{0,1,0,0,0,0}), ...
                            subs(gxz_lin,ae,{0,0,1,0,0,0}), ...
                            subs(gxz_lin,ae,{0,0,0,1,0,0}), ...
                            subs(gxz_lin,ae,{0,0,0,0,1,0}), ...
                            subs(gxz_lin,ae,{0,0,0,0,0,1}) ])

% ahora si reemplazo a por lo que en verdad vale 1/sqrt(3)
Bs_sustitutiva = subs(Bs_sustitutiva, a , sym(1/sqrt(3)))

dt_dx = diff(t,xi)*dxi_dx;
Bb = simplify([ subs(dt_dx,ae,{1,0,0,0,0,0}), ...
                subs(dt_dx,ae,{0,1,0,0,0,0}), ...
                subs(dt_dx,ae,{0,0,1,0,0,0}), ...
                subs(dt_dx,ae,{0,0,0,1,0,0}), ...
                subs(dt_dx,ae,{0,0,0,0,1,0}), ...
                subs(dt_dx,ae,{0,0,0,0,0,1}) ])

%% Integro las matrices con una cuadratura de Gauss-Legendre de orden 2
Kb = int(Bb.'*E*I*Bb*dx_dxi,   xi,-1,1);
Ks = int(Bs_sustitutiva.'*G*Aast*Bs_sustitutiva*dx_dxi,xi,-1,1);

disp('Kb = (E*I/(3*L)) * ');   disp(Kb/(E*I/(3*L)))
disp('Ks = (G*Aast)/(9*L) *'); disp(simplify(Ks/((G*Aast)/(9*L))))         

%% -------------------------------------------------------------------------
%% Calcular correctamente los polinomios de las funciones de forma 1D
function N = calc_N(xp, yp, var)
    % se ve verifican los tamaños de los vectores xp y yp
    nx = length(xp);
    ny = length(yp);
    assert(nx == ny, 'Los vectores xp y yp deben tener el mismo tamaño');

    % se calculan los coeficientes de los polinomios
    c = polyfit(xp, yp, nx-1);
    
    % se eliminan los errores en la aproximación numérica, haciendo los
    % coeficientes demasiado pequeños igual a cero
    c(abs(c) < 1e-10) = 0;
    
    % con los coeficientes corregidos se calculan las funciones de forma
    N = poly2sym(c, var);
end