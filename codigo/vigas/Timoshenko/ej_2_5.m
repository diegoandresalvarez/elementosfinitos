% Ejemplo 4.9 Onate (1995)
% Ejemplo 2.5 Onate (2013)

%% Obtener a partir del elemento de viga de Timoshenko con 
% w = pol grado 3
% t = pol grado 2
% obtenga el elemento de viga de EB, imponiendo que gxz = 0

clear, clc

syms xi w1 w2 w3 w4 t1 t4 t5 L
dx_dxi = L/2;
dxi_dx = 2/L;

% -1           -1/3               0               1/3            1 
%  x-------------x----------------x----------------x-------------x--> xi
%  w1            w2                                w3            w4
%  t1                             t5                             t4

%% Funciones de forma Lagrangianas
N1_3 = calc_N([-1 0 1], [1 0 0], xi); % = xi*(xi-1)/2
N5_3 = calc_N([-1 0 1], [0 1 0], xi); % = (1+xi)*(1-xi)
N4_3 = calc_N([-1 0 1], [0 0 1], xi); % = xi*(xi+1)/2

N1_4 = calc_N([-1 -1/3 1/3 1],[1 0 0 0], xi);
N2_4 = calc_N([-1 -1/3 1/3 1],[0 1 0 0], xi);
N3_4 = calc_N([-1 -1/3 1/3 1],[0 0 1 0], xi);
N4_4 = calc_N([-1 -1/3 1/3 1],[0 0 0 1], xi);

%% Se definen w y t
w = N1_4*w1 + N2_4*w2 + N3_4*w3 + N4_4*w4;
t = N1_3*t1 + N5_3*t5 + N4_3*t4;

%% Se calcula gxz = A + B*xi + C*xi^2
gxz = expand(diff(w,xi)*dxi_dx - t);
A = feval(symengine, 'coeff', gxz, xi, 0); % Aquí se esta llamando a la 
B = feval(symengine, 'coeff', gxz, xi, 1); % función "coeff" del MUPAD
C = feval(symengine, 'coeff', gxz, xi, 2);

%% Se hacen A, B y C iguales a cero y se despejan w2, w3 y t5
% al hacer gxz = 0, obtenemos el EF de Euler-Bernoulli
sol = solve(A==0,B==0,C==0, w2,w3,t5);
disp('w2 = '); disp(sol.w2)
disp('w3 = '); disp(sol.w3)
disp('t5 = '); disp(sol.t5)

%% Se reescribe w en funcion de w2 y w3
w = N1_4*w1 + N2_4*sol.w2 + N3_4*sol.w3 + N4_4*w4;

%% Se reescribe t en funcion de t5
t = N1_3*t1 + N5_3*sol.t5 + N4_3*t4;

%% En este caso como gxz=0, entonces t=diff(w,x)
t_forma2 = diff(w,xi)*dxi_dx;

% Observe que ambos procedimientos dieron lo mismo
expand(t - t_forma2)  % da = 0

%% Se define el vector de movimientos nodales del elemento
ae = {w1,t1,w4,t4};

%% Se calculan las matrices Nw y Nt
% Observe que esta es la matriz N de la viga de EB
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
disp('... ambas son diferentes de 1, por lo que la condición de sólido rígido "sum N_i = 1" no se cumple')

%% Se recalcula dt/dx y se calcula la matriz Bb
% Observe que esta es la matriz Bb de la viga de EB
dt_dx = simplify(diff(t,xi)*dxi_dx);
Bb = simplify([ subs(dt_dx,ae,{1,0,0,0}), ...
                subs(dt_dx,ae,{0,1,0,0}), ...
                subs(dt_dx,ae,{0,0,1,0}), ...
                subs(dt_dx,ae,{0,0,0,1}) ])

%% Se recalcula gxz y se calcula la matriz Bs
gxz = simplify(diff(w,xi)*dxi_dx - t);
%{
Bs = simplify([ subs(gxz,ae,{1,0,0,0,0}), ...
                subs(gxz,ae,{0,1,0,0,0}), ...
                subs(gxz,ae,{0,0,1,0,0}), ...
                subs(gxz,ae,{0,0,0,1,0}), ...
                subs(gxz,ae,{0,0,0,0,1}) ])
%}
% Como gxz == 0, entonces 
Bs = sym([ 0 0 0 0 0 ])

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
    
    % se eliminan los errores en la aproximación numérica, haciendo los
    % coeficientes demasiado pequeños igual a cero
    c(abs(c) < 1e-10) = 0;
    
    % con los coeficientes corregidos se calculan las funciones de forma
    N = poly2sym(c, var);
end
