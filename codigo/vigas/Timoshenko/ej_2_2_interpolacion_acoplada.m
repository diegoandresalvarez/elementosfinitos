% Ejemplo 2.3 Oñate (2013)

%% Interpolación acoplada
clear, clc

syms xi w1 w2 w3 t1 t2 t3 L

dx_dxi = L/2;
dxi_dx = 2/L;

%% Funciones de forma Lagrangianas
N1 = calc_N([-1 0 1], [1 0 0], xi); % = xi*(xi-1)/2
N2 = calc_N([-1 0 1], [0 1 0], xi); % = (1+xi)*(1-xi)
N3 = calc_N([-1 0 1], [0 0 1], xi); % = xi*(xi+1)/2

%% Se definen w y t (w3, t3 es el nodo central)
% la idea de esta numeración es que eliminaremos el nodo central y lo
% expresaremos en función de los gdl asociados a los nodos de los lados
w = N1*w1 + N2*w3 + N3*w2;  % OJO al orden de los términos!!!
t = N1*t1 + N2*t3 + N3*t2;  % OJO al orden de los términos!!!

%% Se calcula gxz = A + B*xi + C*xi^2
gxz = simplify(diff(w,xi)*dxi_dx - t)
A = feval(symengine, 'coeff', gxz, xi, 0); % Aquí se esta llamando a la 
B = feval(symengine, 'coeff', gxz, xi, 1); % función "coeff" del MUPAD
C = feval(symengine, 'coeff', gxz, xi, 2);
% Nota: también se pudo haber usado:
%{
tmp = coeffs(gxz, xi);
A = tmp(1);
B = tmp(2);
C = tmp(3);
%}

%% Se hacen B y C iguales a cero y se despejan w3 y t3, dejando gxz = A = const
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
    
    % se eliminan los errores en la aproximación numérica, haciendo los
    % coeficientes demasiado pequeños igual a cero
    c(abs(c) < 1e-10) = 0;
    
    % con los coeficientes corregidos se calculan las funciones de forma
    N = poly2sym(c, var);
end