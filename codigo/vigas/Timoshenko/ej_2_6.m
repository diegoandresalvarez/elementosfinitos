% Ejemplo 4.10 Oñate (1995)
% Ejemplo 2.6 Oñate (2013)

%% A partir del elemento de viga de Timoshenko con 
% w = pol grado 2
% t = pol grado 2
% e imponiendo que gxz(+/- 1/sqrt(3)) = 0
% obtenga un elemento de viga de dos nodos
% Verifique que la matriz de rigidez del nuevo elemento coincide con la del
% elemento de dos nodos de Euler-Bernoulli

clear, clc

syms xi w1 w2 w3 t1 t2 t3 E I Aast G L
dx_dxi = L/2;
dxi_dx = 2/L;

% -1                              0                              1 
%  x------------------------------x------------------------------x--> xi
%  w1                             w2                             w3
%  t1                             t2                             t3

%% Funciones de forma Lagrangianas
N1_3 = calc_N([-1 0 1], [1 0 0], xi); % = xi*(xi-1)/2
N2_3 = calc_N([-1 0 1], [0 1 0], xi); % = (1+xi)*(1-xi)
N3_3 = calc_N([-1 0 1], [0 0 1], xi); % = xi*(xi+1)/2

%% Se definen w y t
w = N1_3*w1 + N2_3*w2 + N3_3*w3;
t = N1_3*t1 + N2_3*t2 + N3_3*t3;

%% Se calcula gxz
gxz = expand(diff(w,xi)*dxi_dx - t);

%% Se calcula w2 y t2 de modo que gxz(+/- 1/sqrt(3)) = 0
sol = solve(subs(gxz,xi,-1/sqrt(3)) == 0, ...
            subs(gxz,xi,+1/sqrt(3)) == 0, ...   
            w2, t2);

disp('w2 = '); disp(simplify(sol.w2))
disp('t2 = '); disp(simplify(sol.t2))

%% Se reemplaza w2 en la definicion original de w
w = N1_3*w1 + N2_3*sol.w2 + N3_3*w3;
w = collect(w,{'w1', 't1', 'w3', 't3'});

%% Se reemplaza t2 en la definicion original de t
t = N1_3*t1 + N2_3*sol.t2 + N3_3*t3;
t = collect(t,{'w1', 't1', 'w3', 't3'});

%% Se define el vector de movimientos nodales del elemento
ae = {w1,t1,w3,t3};

%% Se calculan las matrices Nw y Nt
Nw = simplify([ subs(w,ae,{1,0,0,0}), ...
                subs(w,ae,{0,1,0,0}), ...
                subs(w,ae,{0,0,1,0}), ...
                subs(w,ae,{0,0,0,1}) ])

Nt = simplify([ subs(t,ae,{1,0,0,0}), ...
                subs(t,ae,{0,1,0,0}), ...
                subs(t,ae,{0,0,1,0}), ...
                subs(t,ae,{0,0,0,1}) ])

%% Se verifica la condicion de cuerpo rigido (sum N_i = 1)
fprintf('sum(Nw) = %s\n', char(expand(sum(Nw))));
% NOTA: en Oñate dice que no cumple con esta condición, sin embargo aquí vemos
% que si cumple, por lo que el libro de Oñate, como cosa rara, tiene un error

fprintf('sum(Nt) = %s\n', char(expand(sum(Nt))));
% Esta condición no se cumple

%% Se recalcula dt/dx y se calcula la matriz Bb
% Observe que esta es la matriz Bb de la viga de EB
dt_dx = simplify(diff(t,xi)*dxi_dx);
Bb = simplify([ subs(dt_dx,ae,{1,0,0,0}), ...
                subs(dt_dx,ae,{0,1,0,0}), ...
                subs(dt_dx,ae,{0,0,1,0}), ...
                subs(dt_dx,ae,{0,0,0,1}) ])
disp('Observe la variacion lineal de Bb (y por lo tanto del momento flector)')

%% Se recalcula gxz y se calcula la matriz Bs
gxz = simplify(diff(w,xi)*dxi_dx - t);
Bs = simplify([ subs(gxz,ae,{1,0,0,0}), ...
                subs(gxz,ae,{0,1,0,0}), ...
                subs(gxz,ae,{0,0,1,0}), ...
                subs(gxz,ae,{0,0,0,1}) ])

% Observe que las raices de gxz son +/- 1/sqrt(3)
assume(L > 0); % esto es para evitar un warning de la siguiente línea
solve(gxz == 0,xi)

%% Se calculan la matriz de rigidez Kb
% esta es la matriz de rigidez K por la teoría de EB
Kb = int(Bb.'*E*I*Bb*L/2, xi,-1,1)

%% Se calculan la matriz de rigidez Ks
% 1) La exacta
Ks_exacta = int(Bs.'*G*Aast*Bs*L/2,xi,-1,1) % La exacta

% 2) Con una cuadratura de GL de orden 2
w = 1; % para xi = -sqrt(1/3), sqrt(1/3)
Ks_GL2 = simplify(subs(Bs.'*G*Aast*Bs*L/2,xi,-sqrt(1/3))*w + ...
                  subs(Bs.'*G*Aast*Bs*L/2,xi,+sqrt(1/3))*w)

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