% Ejemplo 4.8 Onate (1995)
% Ejemplo 2.4 Onate (2013)

% Obtener a partir del elemento de viga de Timoshenko con 
% w = pol grado 2
% t = pol grado 2
% un elemento de variacion lineal para el momento flector y 
% constante para el esfuerzo cortante

% -1                              0                              1 
%  x------------------------------x------------------------------x--> xi
%  w1                             w2                             w3
%  t1                             t2                             t3

clear, clc

syms xi w1 w2 w3 t1 t2 t3 L

dxi_dx = 2/L;

%% Funciones de forma Lagrangianas
N1_2 = poly2sym(polyfit([-1 0 1],[1 0 0],2),xi);   % = xi*(xi-1)/2;
N2_2 = poly2sym(polyfit([-1 0 1],[0 1 0],2),xi);   % = (1+xi)*(1-xi);
N3_2 = poly2sym(polyfit([-1 0 1],[0 0 1],2),xi);   % = xi*(xi+1)/2;

%% Se definen w y t
w = N1_2*w1 + N2_2*w2 + N3_2*w3;
t = N1_2*t1 + N2_2*t2 + N3_2*t3;

%% Se calcula gxz = A + B*xi + C*xi^2
gxz = expand(diff(w,xi)*dxi_dx - t);

A = feval(symengine, 'coeff', gxz, xi, 0); % Aqui se esta llamando a la 
B = feval(symengine, 'coeff', gxz, xi, 1); % funcion "coeff" del MUPAD
C = feval(symengine, 'coeff', gxz, xi, 2);

%% Se hacen B y C iguales a cero y se despejan w2 y t2
% Esto con el animo de usar un campo de deformaciones constante
sol = solve(B==0,C==0, w2,t2);
disp('w2 = '); disp(sol.w2)
disp('t2 = '); disp(sol.t2)

%% Se recalcula dt_dx
tt = N1_2*t1 + N2_2*sol.t2 + N3_2*t3;
dt_dx = simplify(diff(tt,xi)*dxi_dx)
disp('En este caso observe que t2 es una condicion que no nos da un dt_dx')
disp('lineal, sino constante: esto no se pidio en el enunciado')

clear sol tt

%% Se hacen B igual a cero y se despeja w2
% Esto con el animo de usar un campo de deformaciones cuadratico, pero sin
% terminos lineales
sol.w2 = solve(B==0, w2);     % de modo que gxz = A + C*xi^2
disp('w2 = '); disp(sol.w2)

%% En funcion de w2 se reescribe w
w = N1_2*w1 + N2_2*sol.w2 + N3_2*w3;

%% Se define el vector a
ae = {w1,t1,w3,t3,t2};

%% Se calcula la matriz N
Nw = simplify([ ...
subs(w,ae,{1,0,0,0,0}), ...
subs(w,ae,{0,1,0,0,0}), ...
subs(w,ae,{0,0,1,0,0}), ...
subs(w,ae,{0,0,0,1,0}), ...
subs(w,ae,{0,0,0,0,1}) ])

Nt = simplify([ ...
subs(t,ae,{1,0,0,0,0}), ...
subs(t,ae,{0,1,0,0,0}), ...
subs(t,ae,{0,0,1,0,0}), ...
subs(t,ae,{0,0,0,1,0}), ...
subs(t,ae,{0,0,0,0,1}) ])

%% Se verifica la condicion de cuerpo rigido (sum N_i = 1)
fprintf('sum(Nw) = %s\n', char(expand(sum(Nw))));
fprintf('sum(Nt) = %s\n', char(expand(sum(Nt))));

%% Se recalcula dt/dx y se calcula la matriz Bb
dt_dx = simplify(diff(t,xi)*dxi_dx);
Bb = simplify([ ...
subs(dt_dx,ae,{1,0,0,0,0}), ...
subs(dt_dx,ae,{0,1,0,0,0}), ...
subs(dt_dx,ae,{0,0,1,0,0}), ...
subs(dt_dx,ae,{0,0,0,1,0}), ...
subs(dt_dx,ae,{0,0,0,0,1}) ])
disp('Observe la variacion lineal de Bb (y por lo tanto del momento flector)')

%% Se recalcula gxz y se calcula la matriz Bs
gxz = simplify(diff(w,xi)*dxi_dx - t);
Bs = simplify([ ...
subs(gxz,ae,{1,0,0,0,0}), ...
subs(gxz,ae,{0,1,0,0,0}), ...
subs(gxz,ae,{0,0,1,0,0}), ...
subs(gxz,ae,{0,0,0,1,0}), ...
subs(gxz,ae,{0,0,0,0,1}) ])

disp('Observe que Bs es cuadratica, pero esta solo se evalua en los puntos de Gauss')

%% Se calcula gxz en los puntos de Gauss
subs(gxz, xi, +1/sqrt(3))
subs(gxz, xi, -1/sqrt(3))

disp('En los puntos de Gauss Bs es constante y en ambos puntos gxz vale lo mismo')
