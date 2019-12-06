function c1_dibujar_barra_deformada_portico_3d(material, V1, V2, carga, qe_loc, ae_loc, matT, esc)
% Esta funcion dibuja el elemento de portico deformado junto con sus 
% respectivos diagramas de fuerza axial, fuerza cortante y momento flector.
%
% El diagrama de momento flector se grafica en el lado opuesto de la fibra
% a tracción
%
% PARAMETROS DE ENTRADA (junto con algunos ejemplos):
% A = area
% E = E
% I = Ix local
% (x1,y1) y (x2,y2) son las coordenadas de los puntos iniciales
% qxloc = @(x) x.^2; % carga en la dir. del eje x local (function handle)
% qyloc = @(x) 0;    % carga en la dir. del eje y local (function handle)
% qe = [ 0.01        % U1, V1, M1 reacciones del nodo 1 en coord. locales
%       -0.01
%        0.04
%       -0.01        % U2, V2, M2 reacciones del nodo 2 en coord. locales
%        0.02
%       -0.07 ];
% ae = [ 0.01        % u1, v1, t1 desplazamientos nodo 1 en coord. locales
%       -0.01
%        0.04
%       -0.01        % u2, v2, t2 desplazamientos nodo 2 en coord. locales
%        0.02
%       -0.07 ];
% esc_def    = 10;  % escalamiento de la deformada
% esc_faxial = 10;  % escalamiento del diagrama de axiales
% esc_V      = 10;  % escalamiento del diagrama de cortantes
% esc_M      = 10;  % escalamiento del diagrama de momentos

%% se definen algunas constantes
A  = material.A; 
E  = material.E;
Iy = material.Iy;
Iz = material.Iz;

qxloc = carga.qxloc; %qxloc{e};  % axial
qyloc = carga.qyloc; %qyloc{e};  % en dir y
qzloc = carga.qzloc; %qzloc{e};  % en dir z
qtloc = carga.qtloc; %qtloc{e};  % torsion

qe = qe_loc; 
ae = ae_loc;

X  = 1; Y  = 2; Z  = 3;

%% resolver la ecuacion diferencial
npuntos = 1001;
L     = norm(V2-V1);
xinit = linspace(0, L, npuntos);
sol   = bvpinit(xinit, zeros(10,1));
sol   = bvp5c(@ecuacion_diferencial, @condiciones_de_apoyo, sol);

%% Calculos intermedios
s       = sol.x;
Vy      = sol.y( 4,:);         % Fuerza cortante [kN]
Mz      = sol.y( 3,:);         % Momento flector [kN/m]
%thetaz = atan(sol.y(2,:));    % Angulo de giro  [rad]
v       = sol.y( 1,:);         % Desplazamiento en Y' de la viga [m]

Vz      = sol.y( 8,:);         % Fuerza cortante [kN]
My      = sol.y( 7,:);         % Momento flector [kN/m]
%thetay = atan(sol.y(6,:));    % Angulo de giro  [rad]
w       = sol.y( 5,:);         % Desplazamiento en Z' de la viga [m]

axial   = sol.y(10,:);         % Fuerza axial [kN]
u       = sol.y( 9,:);         % Desplazamiento en X' de la viga [m]

% T     = sol.y(11,:);         % Momento torsor [kN/m]

%% Dibujar de deformada
figure(2)
pos = matT*[ s + esc.def*u; esc.def*v; esc.def*w ];

xx = pos(X,:) + V1(X);
yy = pos(Y,:) + V1(Y);
zz = pos(Z,:) + V1(Z);

plot3([V1(X) V2(X)], [V1(Y) V2(Y)], [V1(Z) V2(Z)], 'b-', xx, yy, zz, 'r-','LineWidth',2);

%{
%% Dibujar los diagramas de fuerza axial 
figure(3)
pos = matT*[ s; esc.faxial*axial ]; % escalamiento del diagrama

ss = pos(X,:) + V1(X);
aa = pos(Y,:) + V1(Y);

plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 aa y2], 'r-','LineWidth',2);
text(ss(1),   aa(1),   num2str(-qe(X1)));
text(ss(end), aa(end), num2str(+qe(X2)));

%% Dibujar los diagramas de fuerza cortante
figure(4)
pos = matT*[ s; esc.V*V ]; % escalamiento del diagrama

ss = pos(X,:) + x1;
vv = pos(Y,:) + y1;

plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 vv y2], 'r-','LineWidth',2);
text(ss(1),   vv(1),   num2str(+qe(Y1)));
text(ss(end), vv(end), num2str(-qe(Y2)));

%% Dibujar los diagramas de momento flector
figure(5)
pos = matT*[ s; esc.M*M ]; % escalamiento del diagrama

ss = pos(X,:) + x1;
mm = pos(Y,:) + y1;

plot([V1(X) V2(X)], [V1(Y) V2(Y)], 'b-', [x1 ss x2], [y1 mm y2], 'r-','LineWidth',2);
text(ss(1),   mm(1),   num2str(-qe(M1)));
text(ss(end), mm(end), num2str(+qe(M2)));
[minM,idminM] = min(M); text(ss(idminM), mm(idminM), num2str(minM));
[maxM,idmaxM] = max(M); text(ss(idmaxM), mm(idmaxM), num2str(maxM));
%}

%% ------------------------------------------------------------------------
   function dydx = ecuacion_diferencial(x,y)
      % aqui se implementa la ecuacion diferencial para vigas de material
      % homogeneo y seccion transversal constante (A, E, Iy, Iz, qx, qy, qz
      % las provee la funcion exterior)

      dydx = zeros(10,1);
      %         y(1)           = v
      dydx(1) = y(2);        % = thetaz        %      d^4 v(x)
      dydx(2) = y(3)/(E*Iz); % = Mz/(EIz)      % EIz ---------- = qyloc(x)
      dydx(3) = y(4);        % = Vy            %        dx^4
      dydx(4) = qyloc(x);    % = qyloc

      %         y(5)           = w
      dydx(5) = y(6);        % = thetay        %      d^4 w(x)
      dydx(6) = y(7)/(E*Iy); % = My/(EIy)      % EIy ---------- = qzloc(x)
      dydx(7) = y(8);        % = Vz            %        dx^4 
      dydx(8) = qzloc(x);    % = qzloc
      
      dydx(9) = y(10)/(A*E); % = u             %      d^2 u(x)
      dydx(10) = -qxloc(x);  % = faxial        % A E ---------- = -qxloc(x)
                                               %        dx^2
   end

%% ------------------------------------------------------------------------

   function res = condiciones_de_apoyo(YL,YR)
      % condiciones de apoyo (cond. de frontera de la ecuacion diferencial)
      u1  = 1; v1 = 2; w1 = 3; t1y = 5; t1z =  6;
      u2  = 7; v2 = 8; w2 = 9; t2y = 11; t2z = 12;
                 
      v_  = 1; tz_  =  2; Mz_ = 3; Vy_ = 4;
      w_  = 5; ty_  =  6; My_ = 7; Vz_ = 8;
      u_  = 9; fax_ = 10;      
      
      res = [ % YL: apoyo izquierdo (LEFT), YR: apoyo derecho (RIGHT)
              YL(u_)  - ae(u1)         % uloc(0)      = u1
              YL(v_)  - ae(v1)         % vloc(0)      = v1
              YL(w_)  - ae(w1)         % wloc(0)      = w1
              YL(ty_) + ae(t1y)        % thetayloc(0) = t1y
              YL(tz_) - ae(t1z)        % thetazloc(0) = t1z              
              
              YR(u_)  - ae(u2)         % uloc(L)      = u2
              YR(v_)  - ae(v2)         % vloc(L)      = v2
              YR(w_)  - ae(w2)         % wloc(L)      = w2
              YR(ty_) + ae(t2y)        % thetayloc(L) = t2y
              YR(tz_) - ae(t2z) ];     % thetazloc(L) = t2z
          
%             YR(tx_) - ae(t2x)        % thetaxloc(L) = t2x                        
%             YL(tx_) - ae(t1x)        % thetaxloc(0) = t1x                        
   end
end
%% bye, bye!
