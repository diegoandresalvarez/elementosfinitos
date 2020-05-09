% Ejemplo 2.3 Onate (2013)

%% Interpolacion acoplada
clear, clc

syms xi w1 w2 w3 t1 t2 t3 L

dx_dxi = L/2;
dxi_dx = 2/L;

%% Funciones de forma Lagrangianas
N1 = poly2sym(polyfit([-1 0 1],[1 0 0],2),xi);   % = xi*(xi-1)/2;
N2 = poly2sym(polyfit([-1 0 1],[0 1 0],2),xi);   % = (1+xi)*(1-xi);
N3 = poly2sym(polyfit([-1 0 1],[0 0 1],2),xi);   % = xi*(xi+1)/2;

%% Se definen w y t (w3, t3 es el nodo central)
% la idea de esta numeracion es que eliminaremos el nodo central y lo
% expresaremos en funcion de los gdl asociados a los nodos de los lados
w = N1*w1 + N2*w3 + N3*w2;  % OJO al orden de los terminos!!!
t = N1*t1 + N2*t3 + N3*t2;  % OJO al orden de los terminos!!!

%% Se calcula gxz = A + B*xi + C*xi^2
gxz = simplify(diff(w,xi)*dxi_dx - t)
A = feval(symengine, 'coeff', gxz, xi, 0); % Aqui se esta llamando a la 
B = feval(symengine, 'coeff', gxz, xi, 1); % funcion "coeff" del MUPAD
C = feval(symengine, 'coeff', gxz, xi, 2);
% Nota: tambien se pudo haber usado:
%{
gxz = simplify(diff(w,xi)*dxi_dx - t)
tmp = coeffs(gxz, xi);
A = tmp(1);
B = tmp(2);
C = tmp(3);
%}

%% Se hacen B y C iguales a cero y se despejan w3 y t3
sol = solve(B==0, C==0, w3,t3);
disp('w3 = '); disp(sol.w3)
disp('t3 = '); disp(sol.t3)

%% Y en funcion de w3 y t3 se reescriben w y t
w = N1*w1 + N2*sol.w3 + N3*w2;
w = collect(w, {'w1', 'w2'})

t = N1*t1 + N2*sol.t3 + N3*t2;
t = collect(t, {'t1', 't2'})

%% Se define el vector de movimientos nodales del elemento
ae = {w1,t1,w2,t2};

%% Se calculan las matrices Nw y Nt
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
fprintf('sum(Nt) = %s\n', char(expand(sum(Nt))));

%% Se recalcula dt/dx y se calcula la matriz Bb
dt_dx = simplify(diff(t,xi)*dxi_dx);
Bb = simplify([ subs(dt_dx,ae,{1,0,0,0}), ...
                subs(dt_dx,ae,{0,1,0,0}), ...
                subs(dt_dx,ae,{0,0,1,0}), ...
                subs(dt_dx,ae,{0,0,0,1}) ])
disp('Observe que Bb es constante (y el momento flector tambien)')

%% Se recalcula gxz y se calcula la matriz Bs
gxz = simplify(diff(w,xi)*dxi_dx - t);
Bs = simplify([ subs(gxz,ae,{1,0,0,0}), ...
                subs(gxz,ae,{0,1,0,0}), ...
                subs(gxz,ae,{0,0,1,0}), ...
                subs(gxz,ae,{0,0,0,1}) ])
disp('Observe que Bs es constante (y la fuerza cortante tambien)')

%% Se calculan las matrices de rigidez
syms E I G Aast
Kb = int(Bb.'*E*I*Bb*dx_dxi,   xi,-1,1);
Ks = int(Bs.'*G*Aast*Bs*dx_dxi,xi,-1,1);

disp('Kb = (E*I/L) * '),    pretty(Kb/(E*I/L))
disp('Ks = (G*Aast/L) * '), pretty(Ks/(G*Aast/L))

%% bye, bye!
return
