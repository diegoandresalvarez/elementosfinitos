% Obtener a partir del elemento de viga de Timoshenko con 
% w = pol grado 3
% t = pol grado 2
% un elemento de variacion lineal para el momento flector y 
% constante para el esfuerzo cortante

clear, clc

syms xi w1 w2 w3 w4 t1 t4 t5 L

% -1        -1/sqrt(3)            0            1/sqrt(3)         1 
%  x-------------x----------------x----------------x-------------x--> xi
%  w1            w2                                w3            w4
%  t1                             t5                             t4


dxi_dx = 2/L;

%% Funciones de forma Lagrangianas
N1_2 = poly2sym(polyfit([-1 0 1],[1 0 0],2),xi);   % = xi*(xi-1)/2;
N2_2 = poly2sym(polyfit([-1 0 1],[0 1 0],2),xi);   % = (1+xi)*(1-xi);
N3_2 = poly2sym(polyfit([-1 0 1],[0 0 1],2),xi);   % = xi*(xi+1)/2;

N1_3 = poly2sym(polyfit([-1 -1/3 1/3 1],[1 0 0 0], 3),xi);
N2_3 = poly2sym(polyfit([-1 -1/3 1/3 1],[0 1 0 0], 3),xi);
N3_3 = poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 1 0], 3),xi);
N4_3 = poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 0 1], 3),xi);

%% Se definen w y t
w = N1_3*w1 + N2_3*w2 + N3_3*w3 + N4_3*w4;
t = N1_2*t1 + N2_2*t5 + N3_2*t4;

%% Se calcula gxz = A + B*xi + C*xi^2
gxz = expand(diff(w,xi)*dxi_dx - t);

A = feval(symengine, 'coeff', gxz, xi, 0); % Aqui se esta llamando a la 
B = feval(symengine, 'coeff', gxz, xi, 1); % funcion "coeff" del MUPAD
C = feval(symengine, 'coeff', gxz, xi, 2);

%% Se hacen B y C iguales a cero y se despejan w2 y w3
sol = solve(B,C, w2,w3);
disp('w2 = '); disp(sol.w2)
disp('w3 = '); disp(sol.w3)

%% En funcion de w2 y w3 se reescribe w
w = N1_3*w1 + N2_3*sol.w2 + N3_3*sol.w3 + N4_3*w4;

%% Se define el vector de movimientos nodales del elemento
ae = {w1,t1,w4,t4,t5};

%% Se calcula la matriz N
N = simple([ ...
subs(w,ae,{1,0,0,0,0}), ...
subs(w,ae,{0,1,0,0,0}), ...
subs(w,ae,{0,0,1,0,0}), ...
subs(w,ae,{0,0,0,1,0}), ...
subs(w,ae,{0,0,0,0,1}) ])

%% Se verifica la condicion de cuerpo rigido (sum N_i = 1)
fprintf('sum(N) = %s\n', char(expand(sum(N))));

%% Se recalcula dt/dx y se calcula la matriz Bf
dt_dx = simple(diff(t,xi)*dxi_dx);
Bf = simple([ ...
subs(dt_dx,ae,{1,0,0,0,0}), ...
subs(dt_dx,ae,{0,1,0,0,0}), ...
subs(dt_dx,ae,{0,0,1,0,0}), ...
subs(dt_dx,ae,{0,0,0,1,0}), ...
subs(dt_dx,ae,{0,0,0,0,1}) ])
disp('Observe la variacion lineal de Bf (y por lo tanto del momento flector)')

%% Se recalcula gxz y se calcula la matriz Bc
gxz = simple(diff(w,xi)*dxi_dx - t);
Bc = simple([ ...
subs(gxz,ae,{1,0,0,0,0}), ...
subs(gxz,ae,{0,1,0,0,0}), ...
subs(gxz,ae,{0,0,1,0,0}), ...
subs(gxz,ae,{0,0,0,1,0}), ...
subs(gxz,ae,{0,0,0,0,1}) ])
disp('Observe que Bc es constante (y por lo tanto la fuerza cortante lo es)')