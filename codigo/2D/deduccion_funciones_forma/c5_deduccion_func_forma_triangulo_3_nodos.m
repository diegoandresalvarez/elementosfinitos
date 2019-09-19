%% Programa para deducir las funciones de forma del elemento triangular de
%% tres nodos
clear, clc

%% Se definen las variables simbolicas
syms alpha_1 alpha_2 alpha_3 x1 x2 x3 y1 y2 y3 u1 u2 u3 x y Area

%% Resuelve el sistema de ecuaciones por alpha_1, alpha_2 y alpha_3
r = solve('u1 = alpha_1 + alpha_2*x1 + alpha_3*y1',...
          'u2 = alpha_1 + alpha_2*x2 + alpha_3*y2',...
          'u3 = alpha_1 + alpha_2*x3 + alpha_3*y3',...
          'alpha_1','alpha_2','alpha_3');

%% Muestra alpha_1, alpha_2 y alpha_3
disp('alpha_1 = '); pretty(r.alpha_1)
disp('alpha_2 = '); pretty(r.alpha_2)
disp('alpha_3 = '); pretty(r.alpha_3)

%% Establece de nuevo u y se muestra
u = r.alpha_1 + r.alpha_2*x + r.alpha_3*y;
disp('u = '); pretty(u)

%% Se separa u en su numerador y su denominador
[num,den] = numden(u);

%% Nosotros sabemos que el denominador es dos veces el area del triangulo.
% A continuacion se verificara:
A = det([ 1 x1 y1          % Area del triangulo con vertices
          1 x2 y2          % (x1,y1), (x2,y2) y (x3,y3) numerados en el
          1 x3 y3])/2;     % sentido horario de las manecillas del reloj

disp(simple(den == 2*A))   % el TRUE en la respuesta confirma que el determinante es 2*Area

%% Se vuelve a escribir u, pero con el denominador expresado como 2*Area
% ya que en el paso anterior establecimos que den es igual a 2*A
u = subs(u,den,2*Area);

%% Se factoriza u1, u2 y u3
u = collect(u,u1);
u = collect(u,u2);
u = collect(u,u3);

%% Se muestra finalmente u
% Recuerde que 
% u(x,y) = N1(x,y) u1 + N2(x,y) u2 + N3(x,y) u3
% en la siguiente expresion se pueden ver claramente los terminos de N1, N2
% y N3
disp('u = '); pretty(u)
