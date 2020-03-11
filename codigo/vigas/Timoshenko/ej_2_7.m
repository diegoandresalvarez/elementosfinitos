% Ejemplo 4.11 Oñate (1995)
% Ejemplo 2.7 Oñate (2013)

% A partir del elemento de viga de Timoshenko con 
% w = pol grado 2
% t = pol grado 1
% e imponiendo gxz constante y eliminando el nodo w2
% obtenga un elemento de viga de Timoshenko dos nodos

clear, clc

syms xi w1 w2 w3 t1 t2 t3 E I Aast G L

% -1                              0                              1 
%  x------------------------------x------------------------------x--> xi
%  w1                             w2                             w3
%  t1                                                            t3

dxi_dx = 2/L;

%% Funciones de forma Lagrangianas
N1_1 = poly2sym(polyfit([-1 1],[1 0],1),xi);       % = (1-xi)/2;
N2_1 = poly2sym(polyfit([-1 1],[0 1],1),xi);       % = (xi+1)/2;

N1_2 = poly2sym(polyfit([-1 0 1],[1 0 0],2),xi);   % = xi*(xi-1)/2;
N2_2 = poly2sym(polyfit([-1 0 1],[0 1 0],2),xi);   % = (1+xi)*(1-xi);
N3_2 = poly2sym(polyfit([-1 0 1],[0 0 1],2),xi);   % = xi*(xi+1)/2;

%% Se definen w y t
w = N1_2*w1 + N2_2*w2 + N3_2*w3;
t = N1_1*t1 + N2_1*t3;

%% Se calcula gxz = A + B*xi
gxz = expand(diff(w,xi)*dxi_dx - t);

A = feval(symengine, 'coeff', gxz, xi, 0); % Aqui se esta llamando a la 
B = feval(symengine, 'coeff', gxz, xi, 1); % funcion "coeff" del MUPAD

%% Se hace B igual a cero y se despeja w2
sol.w2 = solve(B==0, w2);
disp('w2 = '); disp(sol.w2)

%% En funcion de w2 se reescribe w
w = N1_2*w1 + N2_2*sol.w2 + N3_2*w3;
w = collect(w,{'w1','t1','w3','t3'});

%% Se define el vector a
ae = {w1,t1,w3,t3};

%% Se calcula la matriz N
N = simplify([ ...
subs(w,ae,{1,0,0,0}), ...
subs(w,ae,{0,1,0,0}), ...
subs(w,ae,{0,0,1,0}), ...
subs(w,ae,{0,0,0,1}) ])

%% Se verifica la condicion de cuerpo rigido (sum N_i = 1)
fprintf('sum(N) = %s\n', char(expand(sum(N))));

%% Se recalcula dt/dx y se calcula la matriz Bb
dt_dx = simplify(diff(t,xi)*dxi_dx);
Bb = simplify([ ...
subs(dt_dx,ae,{1,0,0,0}), ...
subs(dt_dx,ae,{0,1,0,0}), ...
subs(dt_dx,ae,{0,0,1,0}), ...
subs(dt_dx,ae,{0,0,0,1}) ])
disp('Observe la variacion lineal de Bf (y por lo tanto del momento flector)')

%% Se recalcula gxz y se calcula la matriz Bs
gxz = simplify(diff(w,xi)*dxi_dx - t);
Bs = simplify([ ...
subs(gxz,ae,{1,0,0,0}), ...
subs(gxz,ae,{0,1,0,0}), ...
subs(gxz,ae,{0,0,1,0}), ...
subs(gxz,ae,{0,0,0,1}) ])

%% Se calculan la matriz de rigidez Kb
%Kb = int(Bb.'*E*I*Bb*L/2,   xi,-1,1)
syms Db
Kb = int(Bb.'*Db*Bb*L/2,   xi,-1,1)

%% Se calculan la matriz de rigidez Ks
% 1) La exacta
%Ks_exacta = int(Bs.'*G*Aast*Bs*L/2,xi,-1,1) % La exacta
syms Ds
Ks_exacta = int(Bs.'*Ds*Bs*L/2,xi,-1,1) % La exacta

Kb+ Ks_exacta