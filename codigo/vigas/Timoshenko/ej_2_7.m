% Ejemplo 4.11 Onate (1995)
% Ejemplo 2.7 Onate (2013)

%% A partir del elemento de viga de Timoshenko con 
% w = pol grado 2
% t = pol grado 1
% e imponiendo gxz constante y eliminando el nodo w2
% obtenga un elemento de viga de Timoshenko dos nodos

clear, clc
syms xi w1 w2 w3 t1 t2 t3 E I Aast G L
dx_dxi = L/2;
dxi_dx = 2/L;

% -1                              0                              1 
%  x------------------------------x------------------------------x--> xi
%  w1                             w2                             w3
%  t1                                                            t3

%% Funciones de forma Lagrangianas
N1_2 = calc_N([-1 1], [1 0], xi); % = (1-xi)/2
N2_2 = calc_N([-1 1], [0 1], xi); % = (xi+1)/2

N1_3 = calc_N([-1 0 1], [1 0 0], xi); % = xi*(xi-1)/2
N2_3 = calc_N([-1 0 1], [0 1 0], xi); % = (1+xi)*(1-xi)
N3_3 = calc_N([-1 0 1], [0 0 1], xi); % = xi*(xi+1)/2

%% Se definen w y t
w = N1_3*w1 + N2_3*w2 + N3_3*w3;
t = N1_2*t1 + N2_2*t3;

%% Se calcula gxz = A + B*xi
gxz = expand(diff(w,xi)*dxi_dx - t);

A = feval(symengine, 'coeff', gxz, xi, 0); % Aquí se esta llamando a la 
B = feval(symengine, 'coeff', gxz, xi, 1); % función "coeff" del MUPAD

%% Se hace B igual a cero y se despeja w2, imponiendo así gxz = A = const
sol.w2 = solve(B==0, w2);
disp('w2 = '); disp(sol.w2)

%% Se reescribe w en funcion de w2 
w = N1_3*w1 + N2_3*sol.w2 + N3_3*w3;
w = collect(w,{'w1','t1','w3','t3'});

%% Se define el vector de movimientos nodales del elemento ae
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

%% Se verifica la condición de cuerpo rigido (sum N_i = 1)
fprintf('sum(Nw) = %s\n', char(expand(sum(Nw))));
fprintf('sum(Nt) = %s\n', char(expand(sum(Nt))));

%% Se recalcula dt/dx y se calcula la matriz Bb
dt_dx = simplify(diff(t,xi)*dxi_dx);
Bb = simplify([ subs(dt_dx,ae,{1,0,0,0}), ...
                subs(dt_dx,ae,{0,1,0,0}), ...
                subs(dt_dx,ae,{0,0,1,0}), ...
                subs(dt_dx,ae,{0,0,0,1}) ])
disp('Observe que Bb es constante (y el momento flector también)')

%% Se recalcula gxz y se calcula la matriz Bs
gxz = simplify(diff(w,xi)*dxi_dx - t);
Bs = simplify([ subs(gxz,ae,{1,0,0,0}), ...
                subs(gxz,ae,{0,1,0,0}), ...
                subs(gxz,ae,{0,0,1,0}), ...
                subs(gxz,ae,{0,0,0,1}) ])
disp('Observe que Bs es constante (y la fuerza cortante también)')          

%% Se calculan las matrices de rigidez
syms E I G Aast
Kb = int(Bb.'*   E*I*Bb*dx_dxi, xi,-1,1);
Ks = int(Bs.'*G*Aast*Bs*dx_dxi, xi,-1,1);

disp('Kb = (E*I/L) * '),    pretty(Kb/(E*I/L))
disp('Ks = (G*Aast/L) * '), pretty(Ks/(G*Aast/L))

%% bye, bye!
return

%% -------------------------------------------------------------------------
%% Calcular correctamente los polinomios de las funciones de forma 1D
function N = calc_N(xp, yp, var)
    % se ve verifican los tamaños de los vectores xp y yp
    nx = length(xp);
    ny = length(yp);
    assert(nx == ny, 'Los vectores xp y yp deben tener el mismo tamaño');

    % se calculan los coeficientes de los polinomios
    c = polyfit(xp, yp, nx-1);
    
    % se eliminan los errores en la aproximación numerica, haciendo los
    % coeficientes demasiado pequeños igual a cero
    c(abs(c) < 1e-10) = 0;
    
    % con los coeficientes corregidos se calculan las funciones de forma
    N = poly2sym(c, var);
end
