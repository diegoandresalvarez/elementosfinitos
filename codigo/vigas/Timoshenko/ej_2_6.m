% Ejemplo 4.10 Onate (1995)
% Ejemplo 2.6 Onate (2013)

% A partir del elemento de viga de Timoshenko con 
% w = pol grado 2
% t = pol grado 2
% e imponiendo que gxz(+/- 1/sqrt(3)) = 0
% obtenga un elemento de viga de dos nodos
% Verifique que la matriz de rigidez del nuevo elemento coincide con la del
% elemento de dos nodos de Euler-Bernoulli

clear, clc

syms xi w1 w2 w3 t1 t2 t3 E I Aast G L

% -1                              0                              1 
%  x------------------------------x------------------------------x--> xi
%  w1                             w2                             w3
%  t1                             t2                             t3

dxi_dx = 2/L;

%% Funciones de forma Lagrangianas
N1_2 = poly2sym(polyfit([-1 0 1],[1 0 0],2),xi);   % = xi*(xi-1)/2;
N2_2 = poly2sym(polyfit([-1 0 1],[0 1 0],2),xi);   % = (1+xi)*(1-xi);
N3_2 = poly2sym(polyfit([-1 0 1],[0 0 1],2),xi);   % = xi*(xi+1)/2;

%% Se definen w y t
w = N1_2*w1 + N2_2*w2 + N3_2*w3;
t = N1_2*t1 + N2_2*t2 + N3_2*t3;

%% Se calcula gxz
gxz = expand(diff(w,xi)*dxi_dx - t);

%% Se calcula w2 y t2 de modo que gxz(+/- 1/sqrt(3)) = 0
sol = solve(...
   subs(gxz,xi,-1/sqrt(3)) == 0, ...
   subs(gxz,xi,+1/sqrt(3)) == 0, ...   
   w2, t2);

disp('w2 = '); disp(simplify(sol.w2))
disp('t2 = '); disp(simplify(sol.t2))

%% En funcion de w2 se reescribe w
w = N1_2*w1 + N2_2*sol.w2 + N3_2*w3;
w = collect(w,{'w1', 't1', 'w3', 't3'});

t = N1_2*t1 + N2_2*sol.t2 + N3_2*t3;
t = collect(t,{'w1', 't1', 'w3', 't3'});

%% Se define el vector a
ae = {w1,t1,w3,t3};

%% Se calcula la matriz N
Nw = simplify([ ...
subs(w,ae,{1,0,0,0}), ...
subs(w,ae,{0,1,0,0}), ...
subs(w,ae,{0,0,1,0}), ...
subs(w,ae,{0,0,0,1}) ])

Nt = simplify([ ...
subs(t,ae,{1,0,0,0}), ...
subs(t,ae,{0,1,0,0}), ...
subs(t,ae,{0,0,1,0}), ...
subs(t,ae,{0,0,0,1}) ])

%% Se verifica la condicion de cuerpo rigido (sum N_i = 1)
fprintf('sum(Nw) = %s\n', char(expand(sum(Nw))));
% NOTA: en Onate dice que no cumple con esta condicion, sin embargo aqui
% vemos que si cumple, por lo que el libro de Onate tiene un error

fprintf('sum(Nt) = %s\n', char(expand(sum(Nt))));
% Esta condición no se cumple

%% Se recalcula dt/dx y se calcula la matriz Bb
% Observe que esta es la matriz Bb de la viga de EB
dt_dx = simplify(diff(t,xi)*dxi_dx);
Bb = simplify([ ...
subs(dt_dx,ae,{1,0,0,0}), ...
subs(dt_dx,ae,{0,1,0,0}), ...
subs(dt_dx,ae,{0,0,1,0}), ...
subs(dt_dx,ae,{0,0,0,1}) ])
disp('Observe la variacion lineal de Bb (y por lo tanto del momento flector)')

%% Se recalcula gxz y se calcula la matriz Bs
gxz = simplify(diff(w,xi)*dxi_dx - t);
Bs = simplify([ ...
subs(gxz,ae,{1,0,0,0}), ...
subs(gxz,ae,{0,1,0,0}), ...
subs(gxz,ae,{0,0,1,0}), ...
subs(gxz,ae,{0,0,0,1}) ])

% Observe que las raices de gxz son +/- 1/sqrt(3)
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
