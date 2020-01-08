function c4_escamilla_ej_5_5_EB_eq_diff
% Este programa calcula los diagramas de cortante, momento, angulo de giro 
% y deflexion vertical de la viga Ej 5.5 (analisis de estructuras - J. Uribe)
% a partir de la solucion de la ecuacion diferencial:
%
%  d^2  /            d^2v(x)  \
% ------| E(x) I(x) --------- | = q(x)
%  dx^2 \             dx^2    /
%
% utilizando el comando de MATLAB bvp5c

close all; clc;

%% Datos:
b = 0.30;         % Ancho de la viga, m
h = 1.50;         % Altura de la viga, m
E = 200e6;        % Modulo de elasticidad de la viga, GPa
I = (b*h^3)/12;   % Momento de inercia y, m^4

%% resolver la ecuacion diferencial
% Nota: observe que al correr este programa sale el error:
% Warning: Unable to meet the tolerance without using more than 2500 mesh points.
%  The last mesh of 1903 points and the solution are available in the output argument.
%  The maximum error is 0.253843, while requested accuracy is 0.001. 
% > In bvp5c at 341
%  In c4_escamilla_ej_5_5_EB_eq_diff_ver2 at 26
%
% La solución para este es hacer una malla inicial más fina:
% Por ejemplo: 
% xinit = [0:0.0001:10  10:0.0001:16 16:0.0001:19];
% sin embargo con este valor, no converge aun a la exacta, que es la que da
% el MEF.

xinit = [0:0.01:10  10:0.01:16 16:0.01:19];
sol   = bvpinit(xinit, zeros(4,1));
sol   = bvp5c(@f,@bc,sol);

%% Calculos intermedios
V     = sol.y(4,:);          % Fuerza cortante [kN]
M     = sol.y(3,:);          % Momento flector [kN*m]
theta = atan(sol.y(2,:));    % Angulo de giro  [rad]
v     = sol.y(1,:);          % Desplazamiento vertical de la viga [m]
x     = sol.x;

% Se imprimen los resultados
disp('       x       V(x)      M(x)    theta(x)    v(x)');
disp('      [m]      [kN]     [kN*m]   [rad/1000]  [mm]');
disp([x' V' M' 1e3*theta' 1e3*v']);

apoyo1 = 1;
apoyo2 = find(xinit == 10); apoyo2 = apoyo2(1);
apoyo3 = find(xinit == 16); apoyo3 = apoyo3(1);

fprintf('Reaccion en el apoyo 1 = %g kN\n',    V(apoyo1));
fprintf('Momento  en el apoyo 1 = %g kN m\n', -M(apoyo1));
fprintf('Reaccion en el apoyo 2 = %g kN\n',    V(apoyo2+1)-V(apoyo2));   
fprintf('Reaccion en el apoyo 3 = %g kN\n',    V(apoyo3+1)-V(apoyo3));     

%% Hacer los dibujos
z = zeros(1,length(x));
figure(1);
subplot(2,1,1);   
title('Diagramas de angulo de giro y deflexion vertical')
hold on
plot(x, z, '-k', x, v,'-r','LineWidth',2);
ylabel('desplazamiento v(x) [m]');
grid minor
xlabel('eje x [m]')

subplot(2,1,2);   
plot(x, z, '-k', x, theta,'-r','LineWidth',2);
ylabel('angulo de giro theta(x) [rad]');
grid minor

figure(2);
subplot(2,1,1);   
title('Diagramas de cortante y de momento')
hold on
plot(x, z, '-k', x, 1000*M,'-r','LineWidth',2);
ylabel('momento M(x) [N/m]');
grid minor
xlabel('eje x [m]')

subplot(2,1,2);   
plot(x, z, '-k', x, 1000*V,'-r','LineWidth',2);
ylabel('fuerza cortante V(x) [N]');
grid minor

%% -----------------------------------------------------------------------
   function dydx = f(x,y,region)
      % aqui se implementa la ecuacion diferencial para vigas de material
      % homogeneo y seccion transversal constante (E e I las provee la
      % funcion exterior)
      %  d^2  /            d^2v(x)  \
      % ------| E(x) I(x) --------- | = q(x)
      %  dx^2 \             dx^2    /      
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
              YR(tt,1) - YL(tt,2)    % continuidad de theta(x) en x=10
              YR(MM,1) - YL(MM,2)    % continuidad de M(x)     en x=10
              YR(vv,2)               % v(16)    = 0
              % tramo 3:   YL: rodillo, YR: voladizo
              YL(vv,3)               % v(16) = 0
              YR(tt,2) - YL(tt,3)    % continuidad de theta(x) en x=16
              YR(MM,2) - YL(MM,3)    % continuidad de M(x)     en x=16
              YR(VV,3) - 1.5         % V(L) = -1.5 (cortante)
              YR(MM,3)       ];      % M(L) = 0
   end

%% ------------------------------------------------------------------------
   function qq = q(x)
      % carga aplicada a la viga
      if x <= 4.99
         qq = 0;
      elseif x < 5.01
         qq = -3/0.02;       % carga puntual en tramo 1, kN
      elseif x < 10
         qq = 0;
      elseif x < 16
         qq = -1.2;          % carga distribuida en tramo 2, kN/m
      else % if x < 19
         qq = 0;
      end
   end

end
%%End