function c4_viga_tres_apoyos_eq_diff
% Este programa calcula los diagramas de cortante, momento, ángulo de giro 
% y deflexión vertical de la siguiente viga:
%   _________________________________________________
%   |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | -q
%   |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
%   V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V
%   #################################################
%   o                       o                       o
%  /1\                     /2\                     /3\
% /////                   ooooo                   ooooo
%   |----------L/2----------|----------L/2----------|
%
% a partir de la solución de la ecuacion diferencial
%     d4v(x)
% EI ------- = q(x)
%      dx4
%
% utilizando el comando de MATLAB bvp4c()

%% se cierran las figuras y se borra la pantalla
close all; clc;

%% datos:
b = 0.20;       % [m]   Ancho de la viga
h = 0.40;       % [m]   Altura de la viga
E =  2e6;       % [kPa] Módulo de elasticidad de la viga
L =    6;       % [m]   Longitud de la viga
I = (b*h^3)/12; % [m^4] Momento de inercia y

%% resolver la ecuación diferencial
xinit = [0:0.01:L/2  L/2:0.01:L];
sol = bvpinit(xinit, zeros(4,1));
sol = bvp4c(@f, @bc, sol);

%% Cálculos intermedios
V     = sol.y(4,:);       % [kN]   fuerza cortante
M     = sol.y(3,:);       % [kN*m] momento flector
theta = atan(sol.y(2,:)); % [rad]  ángulo de giro
v     = sol.y(1,:);       % [m]    desplazamiento vertical de la viga
x     = sol.x;

% Se imprimen los resultados
disp('       x       V(x)      M(x)    theta(x)    v(x)');
disp('      [m]      [kN]     [kN*m]    [rad]       [m]');
disp([x' V' M' theta' v']);

n = length(xinit);
fprintf('Reacción en el apoyo 1 = %g kN\n', V(1));
fprintf('Reacción en el apoyo 2 = %g kN\n', V(n/2+1)-V(n/2));
fprintf('Reacción en el apoyo 3 = %g kN\n',-V(end));

%% Hacer los dibujos
figure
subplot(4,1,1)
plot(x, V)
ylabel('fuerza cortante V(x) [kN]');
grid minor

subplot(4,1,2)
plot(x,M)
ylabel('momento flector M(x) [kN/m]');
grid minor

subplot(4,1,3)
plot(x,theta)
ylabel('ángulo de giro theta(x) [rad]');
grid minor

subplot(4,1,4)
plot(x,v)
ylabel('desplazamiento v(x) [m]');
grid minor
xlabel('eje x [m]')

%% -----------------------------------------------------------------------
% funciones anidadas -- E e I las provee la función exterior
%

   function dydx = f(x,y,region)
      % aquí se implementa la ecuación diferencial
      %
      %     d4v(x)
      % EI ------- = q(x)      E e I las provee la función exterior
      %      dx4
      
      dydx = zeros(4,1);
      
      %         y(1)          = v
      dydx(1) = y(2);       % = theta
      dydx(2) = y(3)/(E*I); % = M/(EI)
      dydx(3) = y(4);       % = V
      dydx(4) = q(x);       % = q
   end

%% -----------------------------------------------------------------------

   function res = bc(YL,YR)
      vv = 1; tt = 2; MM = 3; VV = 4;
      % condiciones de frontera (externas e internas)
      res = [ YL(vv,1)            % v(0) = 0
              YL(MM,1)            % M(0) = 0
              YR(vv,1)            % v(L/2-) = 0 (tramo 1)
              YL(vv,2)            % v(L/2+) = 0 (tramo 2)
              YR(tt,1) - YL(tt,2) % continuidad de theta(x) en x=L/2
              YR(MM,1) - YL(MM,2) % continuidad de M(x)     en x=L/2
              YR(vv,2)            % v(L) = 0
              YR(MM,2)         ]; % M(L) = 0
   end

%% -----------------------------------------------------------------------

   function qq = q(x)
      % carga aplicada a la viga
      qq = -5; % [kN/m]
   end

end   % c4_viga_tres_apoyos_eq_diff
