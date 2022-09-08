%% Deduccion de las funciones de forma del elemento triangular de tres nodos
clear, clc

%% Se definen las variables simbolicas
syms a1 a2 a3 x1 x2 x3 y1 y2 y3 u1 u2 u3 x y Area

%% Resuelve el sistema de ecuaciones por a1, a2 y a3
r = solve(u1 == a1 + a2*x1 + a3*y1, ...
          u2 == a1 + a2*x2 + a3*y2, ...
          u3 == a1 + a2*x3 + a3*y3,  a1, a2, a3);

%% Muestra a1, a2 y a3
disp('a1 = '); pretty(r.a1)
disp('a2 = '); pretty(r.a2)
disp('a3 = '); pretty(r.a3)

%% Establece de nuevo u y se muestra
u = r.a1 + r.a2*x + r.a3*y;
disp('u = '); pretty(u)

%% Se separa u en su numerador y su denominador
[num,den] = numden(u);

%% Nosotros sabemos que el denominador es dos veces el area del triangulo.
% A continuacion se verificara:
A = det([ 1 x1 y1          % Area del triangulo con vertices
          1 x2 y2          % (x1,y1), (x2,y2) y (x3,y3) numerados en el
          1 x3 y3])/2;     % sentido horario de las manecillas del reloj

disp(simplify(den == 2*A))   % el TRUE en la respuesta confirma que el determinante es 2*Area

%% Se vuelve a escribir u, pero con el denominador expresado como 2*Area
% ya que en el paso anterior establecimos que den es igual a 2*A
u = num/(2*Area);

%% Se factoriza u1, u2 y u3
u = collect(u, [u1, u2, u3]);

%% Se muestra % u(x,y) = N1(x,y) u1 + N2(x,y) u2 + N3(x,y) u3
disp('u = '); pretty(u)

%% Se extraen las funciones de forma
ae = {u1,u2,u3};
N1 = simplify(subs(u,ae,{1,0,0}))
N2 = simplify(subs(u,ae,{0,1,0}))
N3 = simplify(subs(u,ae,{0,0,1}))

%% Se calcula la matriz de masa consistente
N = @(r,s) subs([ N1  0 N2  0 N3  0
                   0 N1  0 N2  0 N3 ], {x,y,Area}, {r,s,A});

NN = @(a,b) N((1-a-b)*x1 + a*x2 + b*x3, (1-a-b)*y1 + a*y2 + b*y3);

syms rho a b t
M = t*int(int(rho*(NN(a,b).') * NN(a,b) * 2*Area, a,0,1-b), b,0,1);
disp('M = rho*t*Area/12 * ')
disp(12*M/(rho*t*Area))
