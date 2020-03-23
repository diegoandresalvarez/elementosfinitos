% Ejemplo 4.12 Onate (1995)
% Ejemplo 2.8 Onate (2013)

% A partir de la interpolacion
% w = pol grado 1
% t = pol grado 1
% e imponiendo gxz constante 
% obtenga un elemento de viga de Timoshenko dos nodos

clear, clc

syms xi w1 w2 t1 t2 E I Aast G L

% -1                                                             1 
%  x-------------------------------------------------------------x--> xi
%  w1                                                            w2
%  t1                                                            t2

dxi_dx = 2/L;

%% Funciones de forma Lagrangianas
N1 = poly2sym(polyfit([-1 1],[1 0],1),xi);       % = (1-xi)/2;
N2 = poly2sym(polyfit([-1 1],[0 1],1),xi);       % = (xi+1)/2;

%% Se definen w y t
w = N1*w1 + N2*w2;
t = N1*t1 + N2*t2;

%% Se calcula gxz = A + B*xi
gxz = expand(diff(w,xi)*dxi_dx - t);

A = feval(symengine, 'coeff', gxz, xi, 0); % Aqui se esta llamando a la 
B = feval(symengine, 'coeff', gxz, xi, 1); % funcion "coeff" del MUPAD

%% Se define el vector a
a = {w1,t1,w2,t2};

%% Se calcula la matriz N
N = simplify([ ...
subs(w,a,{1,0,0,0}), ...
subs(w,a,{0,1,0,0}), ...
subs(w,a,{0,0,1,0}), ...
subs(w,a,{0,0,0,1}) ])

%% Se verifica la condicion de cuerpo rigido (sum N_i = 1)
fprintf('sum(N) = %s\n', char(expand(sum(N))));

%% Se recalcula dt/dx y se calcula la matriz Bb
dt_dx = simplify(diff(t,xi)*dxi_dx);
Bb = simplify([ ...
subs(dt_dx,a,{1,0,0,0}), ...
subs(dt_dx,a,{0,1,0,0}), ...
subs(dt_dx,a,{0,0,1,0}), ...
subs(dt_dx,a,{0,0,0,1}) ])
disp('Observe la variacion lineal de Bb (y por lo tanto del momento flector)')

%% Se recalcula gxz y se calcula la matriz Bs
gxz = simplify(diff(w,xi)*dxi_dx - t);
Bs = simplify([ ...
subs(gxz,a,{1,0,0,0}), ...
subs(gxz,a,{0,1,0,0}), ...
subs(gxz,a,{0,0,1,0}), ...
subs(gxz,a,{0,0,0,1}) ]);

% En vez de hacer B=0, se hace xi=0
Ng = 1;
Bs_sustitutiva = Ng * subs(Bs, xi, 0)

% Se calcula e imprime la matriz de rigidez Ks asociada
Ks = simplify(int(Bs_sustitutiva.'*G*Aast*Bs_sustitutiva*L/2,xi,-1,1));

disp('Ks = (G*Aast/L) *'); disp(simplify(Ks/(G*Aast/L)))