% Ejemplo 4.12 Oñate (1995)
% Ejemplo  2.8 Oñate (2013)

%% A partir de la interpolacion
% w = pol grado 1
% t = pol grado 1
% e imponiendo gxz constante 
% obtenga un elemento de viga de Timoshenko dos nodos

clear, clc

syms xi w1 w2 t1 t2 E I Aast G L
dx_dxi = L/2;
dxi_dx = 2/L;

% -1                                                             1 
%  x-------------------------------------------------------------x--> xi
%  w1                                                            w2
%  t1                                                            t2

%% Funciones de forma Lagrangianas
N1 = calc_N([-1 1], [1 0], xi); % = (1-xi)/2
N2 = calc_N([-1 1], [0 1], xi); % = (1+xi)/2

%% Se definen w y t
w = N1*w1 + N2*w2;
t = N1*t1 + N2*t2;

%% Se calcula gxz = A + B*xi
gxz = expand(diff(w,xi)*dxi_dx - t);
A = feval(symengine, 'coeff', gxz, xi, 0); % Aquí se esta llamando a la 
B = feval(symengine, 'coeff', gxz, xi, 1); % función "coeff" del MUPAD

%% Se define el vector de movimientos nodales ae
ae = {w1,t1,w2,t2};

%% Se calcula la matriz N
Nw = simplify([ subs(w,ae,{1,0,0,0}), ...
                subs(w,ae,{0,1,0,0}), ...
                subs(w,ae,{0,0,1,0}), ...
                subs(w,ae,{0,0,0,1}) ])

Nt = simplify([ subs(t,ae,{1,0,0,0}), ...
                subs(t,ae,{0,1,0,0}), ...
                subs(t,ae,{0,0,1,0}), ...
                subs(t,ae,{0,0,0,1}) ])

%% Se verifica la condición de cuerpo rígido (sum N_i = 1)
fprintf('sum(Nw) = %s\n', char(expand(sum(Nw))));
fprintf('sum(Nt) = %s\n', char(expand(sum(Nt))));

%% Se recalcula dt/dx y se calcula la matriz Bb
dt_dx = simplify(diff(t,xi)*dxi_dx);
Bb = simplify([ subs(dt_dx,ae,{1,0,0,0}), ...
                subs(dt_dx,ae,{0,1,0,0}), ...
                subs(dt_dx,ae,{0,0,1,0}), ...
                subs(dt_dx,ae,{0,0,0,1}) ])
disp('Observe que Bb es constante (y por lo tanto del momento flector)')

%% Se recalcula gxz y se calcula la matriz Bs
gxz = simplify(diff(w,xi)*dxi_dx - t);
Bs = simplify([ subs(gxz,ae,{1,0,0,0}), ...
                subs(gxz,ae,{0,1,0,0}), ...
                subs(gxz,ae,{0,0,1,0}), ...
                subs(gxz,ae,{0,0,0,1}) ]);

% En vez de hacer B=0, se hace xi=0
Ng = 1;
Bs_sustitutiva = Ng * subs(Bs, xi, 0)

% Se calcula e imprime la matriz de rigidez Ks asociada
Ks = simplify(int(Bs_sustitutiva.'*G*Aast*Bs_sustitutiva*L/2,xi,-1,1));

disp('Ks = (G*Aast/L) *'); disp(simplify(Ks/(G*Aast/L)))

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