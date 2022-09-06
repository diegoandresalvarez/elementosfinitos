% Ejemplo 4.8 Onate (1995)
% Ejemplo 2.4 Onate (2013)

%% A partir del elemento de viga de Timoshenko con 
% w = pol grado 2
% t = pol grado 2
% obtener un elemento de variación lineal para el momento flector y 
% constante para el esfuerzo cortante

clear, clc

syms xi w1 w2 w3 t1 t2 t3 L

dx_dxi = L/2;
dxi_dx = 2/L;

% -1                              0                              1 
%  x------------------------------x------------------------------x--> xi
%  w1                             w2                             w3
%  t1                             t2                             t3

%% Funciones de forma Lagrangianas
N1_3 = calc_N([-1 0 1], [1 0 0], xi); % = xi*(xi-1)/2
N2_3 = calc_N([-1 0 1], [0 1 0], xi); % = (1+xi)*(1-xi)
N3_3 = calc_N([-1 0 1], [0 0 1], xi); % = xi*(xi+1)/2

%% Se definen w y t
w = N1_3*w1 + N2_3*w2 + N3_3*w3;
t = N1_3*t1 + N2_3*t2 + N3_3*t3;

%% Se calcula gxz = A + B*xi + C*xi^2
gxz = expand(diff(w,xi)*dxi_dx - t);

A = feval(symengine, 'coeff', gxz, xi, 0); % Aqui se esta llamando a la 
B = feval(symengine, 'coeff', gxz, xi, 1); % función "coeff" del MUPAD
C = feval(symengine, 'coeff', gxz, xi, 2);

%{
% EN EL LIBRO DE ONATE SE EXPLICA QUE ESTO NO FUNCIONA:
%% Se hacen B y C iguales a cero y se despejan w2 y t2
% Esto con el animo de usar un campo de deformaciones constante
sol = solve(B==0,C==0, w2,t2);
disp('w2 = '); disp(sol.w2)
disp('t2 = '); disp(sol.t2)

%% Se recalcula dt_dx
tt = N1_3*t1 + N2_3*sol.t2 + N3_3*t3;
dt_dx = simplify(diff(tt,xi)*dxi_dx)
disp('En este caso observe que t2 es una condicion que no nos da un dt_dx')
disp('lineal, sino constante: esto no se pidio en el enunciado')
%}

%% Se hacen B igual a cero y se despeja w2
% Esto con el animo de usar un campo de deformaciones cuadratico, pero sin
% terminos lineales
sol.w2 = solve(B==0, w2);     % de modo que gxz = A + C*xi^2
disp('w2 = '); disp(sol.w2)

%% En funcion de w2 se reescribe w
w = N1_3*w1 + N2_3*sol.w2 + N3_3*w3;

%% Se define el vector de movimientos nodales ae
ae = {w1,t1,w3,t3,t2};

%% Se calcula las matrices Nw y Nt
Nw = simplify([ subs(w,ae,{1,0,0,0,0}), ...
                subs(w,ae,{0,1,0,0,0}), ...
                subs(w,ae,{0,0,1,0,0}), ...
                subs(w,ae,{0,0,0,1,0}), ...
                subs(w,ae,{0,0,0,0,1}) ])

Nt = simplify([ subs(t,ae,{1,0,0,0,0}), ...
                subs(t,ae,{0,1,0,0,0}), ...
                subs(t,ae,{0,0,1,0,0}), ...
                subs(t,ae,{0,0,0,1,0}), ...
                subs(t,ae,{0,0,0,0,1}) ])

%% Se verifica la condicion de cuerpo rigido (sum N_i = 1)
fprintf('sum(Nw) = %s\n', char(expand(sum(Nw))));
fprintf('sum(Nt) = %s\n', char(expand(sum(Nt))));

%% Se recalcula dt/dx y se calcula la matriz Bb
dt_dx = simplify(diff(t,xi)*dxi_dx);
Bb = simplify([ subs(dt_dx,ae,{1,0,0,0,0}), ...
                subs(dt_dx,ae,{0,1,0,0,0}), ...
                subs(dt_dx,ae,{0,0,1,0,0}), ...
                subs(dt_dx,ae,{0,0,0,1,0}), ...
                subs(dt_dx,ae,{0,0,0,0,1}) ])
disp('Observe la variacion lineal de Bb (y por lo tanto del momento flector)')

%% Se recalcula gxz y se calcula la matriz Bs
gxz = simplify(diff(w,xi)*dxi_dx - t);
Bs = simplify([ subs(gxz,ae,{1,0,0,0,0}), ...
                subs(gxz,ae,{0,1,0,0,0}), ...
                subs(gxz,ae,{0,0,1,0,0}), ...
                subs(gxz,ae,{0,0,0,1,0}), ...
                subs(gxz,ae,{0,0,0,0,1}) ])
disp('Observe que Bs es cuadratica, pero esta solo se evalua en los puntos de Gauss')

%% Se calcula gxz en los puntos de Gauss
subs(gxz, xi, +1/sqrt(3))
subs(gxz, xi, -1/sqrt(3))

disp('En los puntos de Gauss Bs es constante y en ambos puntos gxz vale lo mismo')

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