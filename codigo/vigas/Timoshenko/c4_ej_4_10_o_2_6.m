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
   subs(gxz,xi,-1/sqrt(3)), ...
   subs(gxz,xi,+1/sqrt(3)), ...   
   w2, t2);

disp('w2 = '); disp(sol.w2)
disp('t2 = '); disp(sol.t2)

%% En funcion de w2 se reescribe w
w = N1_2*w1 + N2_2*sol.w2 + N3_2*w3;
w = collect(w,w1);
w = collect(w,t1);
w = collect(w,w3);
w = collect(w,t3)

t = N1_2*t1 + N2_2*sol.t2 + N3_2*t3;
t = collect(t,w1);
t = collect(t,t1);
t = collect(t,w3);
t = collect(t,t3)

%% Se define el vector a
a = {w1,t1,w3,t3};

%% Se calcula la matriz N
N = simplify([ ...
subs(w,a,{1,0,0,0}), ...
subs(w,a,{0,1,0,0}), ...
subs(w,a,{0,0,1,0}), ...
subs(w,a,{0,0,0,1}) ])

%% Se verifica la condicion de cuerpo rigido (sum N_i = 1)
fprintf('sum(N) = %s\n', char(expand(sum(N))));
% NOTA: en Onate dice que no cumple con esta condicion, sin embargo aqui
% vemos que si cumple, por lo que el libro de Onate tiene un error

%% Se recalcula dt/dx y se calcula la matriz Bf
dt_dx = simplify(diff(t,xi)*dxi_dx);
Bf = simplify([ ...
subs(dt_dx,a,{1,0,0,0}), ...
subs(dt_dx,a,{0,1,0,0}), ...
subs(dt_dx,a,{0,0,1,0}), ...
subs(dt_dx,a,{0,0,0,1}) ])
disp('Observe la variacion lineal de Bf (y por lo tanto del momento flector)')

%% Se recalcula gxz y se calcula la matriz Bc
gxz = simplify(diff(w,xi)*dxi_dx - t);
Bc = simplify([ ...
subs(gxz,a,{1,0,0,0}), ...
subs(gxz,a,{0,1,0,0}), ...
subs(gxz,a,{0,0,1,0}), ...
subs(gxz,a,{0,0,0,1}) ])

% Observe que las raices de gxz son +/- 1/sqrt(3)
solve(gxz,xi)

%% Se calculan la matriz de rigidez Kf
Kf = int(Bf.'*E*I*Bf*L/2,   xi,-1,1)

%% Se calculan la matriz de rigidez Kc
% 1) La exacta
Kc_exacta = int(Bc.'*G*Aast*Bc*L/2,xi,-1,1) % La exacta

% 2) Con una cuadratura de GL de orden 2
w = 1; % para xi = -sqrt(1/3), sqrt(1/3)
Kc_GL2 = simplify(subs(Bc.'*G*Aast*Bc*L/2,xi,-sqrt(1/3))*w + subs(Bc.'*G*Aast*Bc*L/2,xi,sqrt(1/3))*w)
