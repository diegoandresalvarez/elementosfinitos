function c4_escamilla_ej_5_5_EB_eq_diff
% Este programa calcula los diagramas de cortante, momento, ángulo de giro 
% y deflexión vertical de la viga Ej 5.5 (análisis de estructuras - J. Uribe)
% a partir de la solución de la ecuacion diferencial:
%
%     d4v(x)
% EI ------- = q(x)      E e I las provee la función exterior
%      dx4
%
% utilizando el comando de MATLAB bvp5c()

%% se cierran las figuras y se borra la pantalla
%close all; clc;

%% datos:
b = 0.30;         % [m]   Ancho de la viga
h = 1.50;         % [m]   Altura de la viga
E = 200e6;        % [kPa] Módulo de elasticidad de la viga
I = (b*h^3)/12;   % [m^4] Momento de inercia y

%% resolver la ecuación diferencial
xinit = [0:0.01:10  10:0.01:16  16:0.01:19];  % [m]
sol   = bvpinit(xinit, zeros(4,1));
opts  = bvpset('NMax', 5000);
sol   = bvp5c(@f, @bc, sol, opts);

% NOTA: opts  = bvpset('NMax', 5000); se utiliza para evitar que salga el
% mensaje:
% Warning: Unable to meet the tolerance without using more than 2500 mesh points.
%  The last mesh of 1903 points and the solution are available in the output argument.
%  The maximum error is 0.253843, while requested accuracy is 0.001. 

%% Cálculos intermedios
V     = sol.y(4,:);       % [kN]   fuerza cortante
M     = sol.y(3,:);       % [kN*m] momento flector
theta = atan(sol.y(2,:)); % [rad]  ángulo de giro
v     = sol.y(1,:);       % [m]    desplazamiento vertical de la viga
x     = sol.x;

% Se imprimen los resultados
disp('       x       V(x)      M(x)    theta(x)    v(x)');
disp('      [m]      [kN]     [kN*m]   [rad/1000]  [mm]');
disp([x' V' M' 1e3*theta' 1e3*v']);

apoyo1 = 1;
apoyo2 = find(x == 10); apoyo2 = apoyo2(1);
apoyo3 = find(x == 16); apoyo3 = apoyo3(1);

fprintf('Reacción en el apoyo 1 = %g kN\n',    V(apoyo1));
fprintf('Momento  en el apoyo 1 = %g kN m\n', -M(apoyo1));
fprintf('Reacción en el apoyo 2 = %g kN\n',    V(apoyo2+1)-V(apoyo2));
fprintf('Reacción en el apoyo 3 = %g kN\n',    V(apoyo3+1)-V(apoyo3));

%% Hacer los dibujos
z = zeros(1, length(x));
figure(1);
subplot(2,1,1);   
title('Diagramas de ángulo de giro y deflexión vertical')
hold on
plot(x, z, '-k', x, v, '-r','LineWidth',2);
ylabel('desplazamiento v(x) [m]');
grid minor
xlabel('eje x [m]')

subplot(2,1,2);
plot(x, z, '-k', x, theta, '-r','LineWidth',2);
ylabel('ángulo de giro theta(x) [rad]');
grid minor
xlabel('eje x [m]')

figure(2);
subplot(2,1,1);
title('Diagramas de fuerza cortante y de momento flector')
hold on
plot(x, z, '-k', x, M,'-r','LineWidth',2);
ylabel('momento flector M(x) [kN/m]');
grid minor
xlabel('eje x [m]')

subplot(2,1,2);
plot(x, z, '-k', x, V,'-r','LineWidth',2);
ylabel('fuerza cortante V(x) [kN]');
grid minor
xlabel('eje x [m]')

%% -----------------------------------------------------------------------
   function dydx = f(x,y,region)
      % aquí se implementa la ecuación diferencial para vigas de material
      % homogeneo y sección transversal constante (E e I las provee la
      % función exterior)
      %
      %     d4v(x)
      % EI ------- = q(x)      E e I las provee la función exterior
      %      dx4
      %
      dydx = zeros(4,1);

      %         y(1)          = v
      dydx(1) = y(2);       % = theta
      dydx(2) = y(3)/(E*I); % = M/(EI)
      dydx(3) = y(4);       % = V
      dydx(4) = q(x);       % = q
   end

%% ------------------------------------------------------------------------
   function res = bc(YL,YR)
      vv = 1; tt = 2; MM = 3; VV = 4;
      % condiciones de frontera (externas e internas)
      res = [ % tramo 1:   YL: empotramiento, YR: rodillo
              YL(vv,1)               % v(0)     = 0
              YL(tt,1)               % theta(0) = 0
              YR(vv,1)               % v(10)    = 0
              % tramo 2:   YL: rodillo, YR: rodillo
              YL(vv,2)               % v(10) = 0
              YR(tt,1) - YL(tt,2)    % continuidad de theta(x)   en x=10 m
              YR(MM,1) - YL(MM,2)    % continuidad de M(x)
              YR(vv,2)               % v(16)    = 0
              % tramo 3:   YL: rodillo, YR: voladizo
              YL(vv,3)               % v(16) = 0
              YR(tt,2) - YL(tt,3)    % continuidad de theta(x)   en x=16 m
              YR(MM,2) - YL(MM,3)    % continuidad de M(x)       
              YR(VV,3) - 15          % V(L) = -15 kN (cortante)  en x=19 m
              YR(MM,3)       ];      % M(L) = 0                  
   end

%% ------------------------------------------------------------------------
   function qq = q(x)
      % carga aplicada a la viga
      if x <= 4.99
         qq = 0;        % x = [ 0.00,  4.99]
      elseif x < 5.01
         qq = -30/0.02; % x = [ 4.99,  5.01]  carga puntual en tramo 1, kN
      elseif x < 10
         qq = 0;        % x = [ 5.01, 10.00]
      elseif x < 16
         qq = -12;      % x = [10.00, 16.00]  carga distribuida en tramo 2, kN/m
      else % if x < 19
         qq = 0;
      end
   end

end   % c4_escamilla_ej_5_5_EB_eq_diff
