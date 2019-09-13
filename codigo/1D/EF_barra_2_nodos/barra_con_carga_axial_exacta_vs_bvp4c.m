%% Definicion del problema
% Calcule los desplazamientos y las reacciones en el empotramiento 
% de la viga mostrada resolviendo la ecuacion diferencial numericamente con
% la funcion bvp4c
% 
% | b (carga distribuida de magnitud b)
% |->->->->->->->->->->->->->->->->
% |================================o--> P (carga puntual P en extremo derecho)
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
% ----| E(x) A(x)------- | + b(x) en x \in [0,L]     dado u(0) = 0
%  dx \            dx    /                                faxial(L) = P
%

% En el caso mas general E, A y b son funciones. Escriba aqui las funciones
% como tal en caso de tener un caso mas general
EE = @(x) E;
AA = @(x) A;
bb = @(x) b;

% Se define la ecuacion diferencial, expresada como un sistema de dos
% ecuaciones diferenciales
u      = 1; % y(1) = u(x)       
faxial = 2; % y(2) = faxial(x)
sist_eq_dif = @(x,y) [ y(faxial)/(EE(x)*AA(x)) 
                       -bb(x)                  ];

% Se definen las condiciones de frontera
% y_izq = condiciones de frontera del lado izquierdo (x=0)
% y_izq(1) = u(x=0)          y_izq(2) = faxial(x=0)
% y_der = condiciones de frontera del lado derecho   (x=L)
% y_der(1) = u(x=L)          y_der(2) = faxial(x=L)
cond_frontera = @ (y_izq,y_der) ...
                 [ y_izq(u)             % u(x=0)      = 0 (desplazamiento)
                   y_der(faxial) - P ]; % faxial(x=L) = P (carga axial)

% Solucion tentativa de la ecuacion diferencial
x = linspace(0,L,30);         % 30 puntos uniformemente distrib. entre 0 y L
y_inicial = bvpinit(x,[0 0]); % el [ 0 0 ] hace y_inicial.y = zeros(2,30)

% Solucion como tal de la ecuacion diferencial
sol = bvp4c(sist_eq_dif, cond_frontera, y_inicial);

% Evaluar la respuesta en los puntos x
y = sol.y;

%% Solucion analitica
u_exacta      = @(x) (-b*x.^2/2 + (P + b*L)*x)/(E*A); % desplazamiento
faxial_exacta = @(x) (P + b*(L-x));                   % carga axial

%% Se reportan los errores en el calculo
error_en_u = max(abs(u_exacta(x) - y(u,:)));
fprintf('Maximo error en el calculo del desplazamiento = %g m\n', error_en_u);

error_en_faxial = max(abs(faxial_exacta(x) - y(faxial,:)));
fprintf('Maximo error en el calculo de la fuerza axial = %g N\n', error_en_faxial);

%% Grafico la solucion analitica y la solucion por el la funcion bvp4c
figure                       % cree un nuevo lienzo

% 1) grafico los desplazamientos de la barra
subplot(2,1,1);              % grafique en la parte superior (1) del lienzo
plot(x, u_exacta(x), 'r');   % grafico solucion analitica
hold on;                     % no borre el lienzo 
plot(x, y(u,:), 'bx');       % grafico solucion por bvp4c()
title('Comparacion de la solucion analitica vs la funcion bvp4c() para el desplazamiento');
xlabel('Eje X (m)')          % titulo del eje X
ylabel('Desplazamiento (m)') % titulo del eje Y
legend('solucion exacta de u(x)','solucion por bvp4c()', ...
                                                   'Location','SouthEast');

% 2) grafico la carga axial de la barra
subplot(2,1,2);              % grafique en la parte inferior (2) del lienzo
plot(x, faxial_exacta(x), 'r');  % grafico solucion analitica
hold on;                     % no borre el lienzo
plot(x, y(faxial,:), 'bx');  % grafico solucion por bvp4c()
title('Comparacion de la solucion analitica vs la funcion bvp4c() para la carga axial');
xlabel('Eje X (m)')          % titulo del eje X
ylabel('Carga axial (N)')    % titulo del eje Y
legend('solucion exacta de f_{axial}(x)','solucion por bvp4c()', ...
                                                   'Location','NorthEast');
%% bye, bye!!!
return;
