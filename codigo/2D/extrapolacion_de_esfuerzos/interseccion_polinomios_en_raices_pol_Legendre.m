%% Programa para ilustrar que en los puntos de una cuadratura de
%% Gauss-Legendre de orden n, un polinomio de grado n y otro de grado n-1
%% obtenido como ajuste por mínimos cuadrados del anterior, toman igual
%% valor.

clear, clc, close all        % borro la memoria, la pantalla y cierro figuras
syms x a b                   % defino las variables simbólicas

% el polinomio f y su aproximación g
f = 1 + x + x^2;  g = a + b*x;  

E = int((f-g)^2,x,-1,1);     % el ajuste por mínimos cuadrados
r = solve(diff(E,a) == 0,...
          diff(E,b) == 0,...
          a, b);
     
g = r.a + r.b*x;             % defino de nuevo la función g
inter = solve(f == g, x)     % calculo la intersección de ambos polinomios
       
figure; hold on; grid on     % grafico la respuesta
fplot([f, g], [-1 1], 'LineWidth', 2);
gg = matlabFunction(g, 'Vars', x);
plot(inter,gg(inter),'r.','MarkerSize',40);
legend('f(x)', 'g(x)', 'intersección','Location','best')
title('Intersección en las raices del polinomio de Legendre de grado 2')

%%
clear                        % borro la memoria
syms x a b c                 % defino las variables simbólicas

% el polinomio f y su aproximación g
f = 1 + x + x^2 + x^3;  g = a + b*x + c*x^2;

E = int((f-g)^2,x,-1,1);     % el ajuste por mínimos cuadrados
r = solve(diff(E,a) == 0,...
          diff(E,b) == 0,...
          diff(E,c) == 0,...
          a, b, c);
     
g = r.a + r.b*x + r.c*x^2;   % defino de nuevo la función g
inter = solve(f == g, x)     % calculo la intersección de ambos polinomios
       
figure; hold on; grid on     % grafico la respuesta
fplot([f, g], [-1 1], 'LineWidth', 2);
gg = matlabFunction(g, 'Vars', x);
plot(inter,gg(inter),'r.','MarkerSize',40);
legend('f(x)', 'g(x)', 'intersección','Location','best')
title('Intersección en las raices del polinomio de Legendre de grado 3')
