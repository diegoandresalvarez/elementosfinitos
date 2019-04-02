clear, clc, close all        % borro la memoria, la pantalla y las figuras

%% Definicion del problema
% Calcule los desplazamientos y las reacciones en el empotramiento 
% de la viga mostrada resolviendo la ecuacion diferencial numericamente con
% la funcion bvp4c
% 
% | b (carga distribuida de magnitud b)
% |->->->->->->->->->->->->->->->->
% |================================o--> P (carga puntual P en nodo nno)
% |<----longitud L de la barra---->|    el area transversal de la barra es A

%% Defino las variables
E = 200e9;    % Pa          % modulo de elasticidad de la barra
A = (0.01)^2; % m^2         % area transversal de la barra
L = 2;        % m           % longitud de la barra
b = 1000;     % N/m         % fuerza axial aplicada sobre cada EF
P = 250;      % N           % carga nodal al final de la barra

%% Solucion de la ecuacion diferencial

% Solucion numerica usando bvp4c (boundary value problem - MATLAB)
%   d /           du(x)  \
% ----| E(x) A(x)------- | + b(x) en x \in [0,L]     dado u(0)=0
%  dx \            dx    /                                faxial(L) = P
%
% bvp4c is a finite difference code that implements the three-stage Lobatto
% IIIa formula. This is a collocation formula and the collocation polynomial
% provides a C1-continuous solution that is fourth-order accurate uniformly 
% in [a,b]. Mesh selection and error control are based on the residual of 
% the continuous solution. 

% En el caso mas general E, A y b son funciones. Escriba aqui las funciones
% como tal en caso de tener un caso mas general
EE = @(x) E;
AA = @(x) A;
bb = @(x) b;

% Por favor, antes de continuar mirar la ayuda de MATLAB de los comandos
% bvpinit, bvp4c, deval

% Se define la ecuacion diferencial, expresada como un sistema de dos
% ecuaciones diferenciales
% y(1) = u(x)
% y(2) = faxial(x)
sist_eq_dif = @(x,y) [ y(2)/(EE(x)*AA(x)) 
                       -bb(x)             ];

% Se definen las condiciones de frontera
% ya = condiciones de frontera del lado izquierdo (x=0)
%     ya(1) = u(x=0)          ya(2) = faxial(x=0)
% yb = condiciones de frontera del lado derecho   (x=L)
%     yb(1) = u(x=L)          yb(2) = faxial(x=L)
cond_frontera = @ (ya,yb) [ ya(1)         % u(x=0)      = 0 (desplazamiento)
                            yb(2) - P ];  % faxial(x=L) = P (carga axial)

% Solucion tentativa de la ecuacion diferencial
x = linspace(0,L,30);           % 30 puntos uniformemente distrib. entre 0 y L
sol_inicial = bvpinit(x,[0 0]); % el [ 0 0 ] hace sol_inicial.y = zeros(2,30)

% Solucion como tal de la ecuacion diferencial
sol = bvp4c(sist_eq_dif, cond_frontera, sol_inicial);

% Evaluar la respuesta en los puntos x
y = deval(sol,x);

%% Solucion analitica
u_exacta      = @(x) (-b*x.^2/2 + (P + b*L)*x)/(E*A); % desplazamiento
faxial_exacta = @(x) (P + b*(L-x));                   % carga axial

%% Grafico la solucion analitica y la solucion por el la funcion bvp4c
figure                             % cree un nuevo lienzo

% 1) grafico los desplazamientos de la barra
subplot(2,1,1);                    % grafique en la parte superior (1) del lienzo
plot(x, u_exacta(x), 'r');         % grafico solucion analitica
hold on;                           % no borre el lienzo 
plot(x, y(1,:), 'bx');             % grafico solucion por bvp4c
title('Comparacion de la solucion analitica vs la funcion bvp4c para el desplazamiento');
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Desplazamiento (m)')       % titulo del eje Y
legend('solucion exacta','solucion por bvp4c', 'Location','SouthEast');

error_en_u = 100*max(abs((u_exacta(x) - y(1,:))./u_exacta(x)));
fprintf('Maximo porcentaje de error en el calculo del desplazamiento = %g%%\n', error_en_u);


% 2) grafico la carga axial de la barra
subplot(2,1,2);                    % grafique en la parte inferior (2) del lienzo
plot(x, faxial_exacta(x), 'r');    % grafico solucion analitica
hold on;                           % no borre el lienzo
plot(x, y(2,:), 'bx');             % grafico solucion por bvp4c
title('Comparacion de la solucion analitica vs la funcion bvp4c para la carga axial');
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Carga axial (N)')          % titulo del eje Y
legend('solucion exacta','solucion por bvp4c', 'Location','NorthEast');

error_en_faxial = 100*max(abs((faxial_exacta(x) - y(2,:))./faxial_exacta(x)));
fprintf('Maximo porcentaje de error en el calculo de la fuerza axial = %g%%\n', error_en_faxial);

%% bye, bye!!!
return;
