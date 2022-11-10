function feglob = calc_fuerzas_nodales_equivalentes(A, E, I, ...
    x1,y1, x2,y2, qxloc,qyloc)
% Este programa calcula las fuerzas nodales equivalentes asociadas a una 
% carga distribuida y retorna los resultados en coordenadas globales.
%
% PARAMETROS:
% A = area
% E = E
% I = Ix local
% (x1,y1) y (x2,y2) son las coordenadas de los puntos iniciales
% qxloc = @(x) x.^2; % carga en la dir. del eje x local (function handle)
% qyloc = @(x) 0;    % carga en la dir. del eje y local (function handle)

% EJEMPLO
%{
A = .30*.35;
E = 190e4;
I = .30*.35^3/12;
x1 = 3; y1 = 4; x2 = 7; y2 = 6;
ang = atan2(y2-y1,x2-x1);
c = cos(ang); s = sin(ang);
qxloc = @(x) -2.8*s*c;
qyloc = @(x) -2.8*c^2;
qeglob = calc_fuerzas_nodales_equivalentes(A, E, I, x1,y1, x2,y2, qxloc,qyloc)
%}

%% se definen algunas constantes
X1  = 1; Y1 = 2; M1 = 3; X2 = 4; Y2 = 5; M2 = 6; 

%% se calcula la longitud de la barra
L = hypot(x2-x1, y2-y1);

%% resolver la ecuacion diferencial
npuntos = 101;
xinit = linspace(0, L, npuntos);
sol   = bvpinit(xinit, zeros(6,1));
sol   = bvp5c(@ecuacion_diferencial, @condiciones_de_apoyo, sol);

%% Calculos intermedios
y = deval(sol, [0 L]);
faxial = y(6,:);           % Fuerza axial [kN]
V      = y(4,:);           % Fuerza cortante [kN]
M      = y(3,:);           % Momento flector [kN/m]
%u     = y(5,:);           % Desplazamiento horizontal de la viga [m]
%v     = y(1,:);           % Desplazamiento vertical de la viga [m]
%theta = atan(sol.y(2,:)); % Angulo de giro  [rad]
X1 = +faxial(1);   Y1 = -V(1);   M1 = +M(1);   % 1   => en x=0
X2 = -faxial(end); Y2 = +V(end); M2 = -M(end); % end => en x=L

feloc = [ X1; Y1; M1; X2; Y2; M2 ];

c = (x2-x1)/L; s = (y2-y1)/L;
T = [ c  s  0  0  0  0        
     -s  c  0  0  0  0        
      0  0  1  0  0  0
      0  0  0  c  s  0
      0  0  0 -s  c  0
      0  0  0  0  0  1];
  
feglob = T'*[ X1; Y1; M1; X2; Y2; M2 ];  

%% -----------------------------------------------------------------------
   function dydx = ecuacion_diferencial(x,y)
      % aqui se implementa la ecuacion diferencial para vigas de material
      % homogeneo y seccion transversal constante (A, E, I, qx, qy las 
      % provee la funcion exterior)
      %      d^4 v(x)
      % E I ---------- = q(x)
      %        dx^4
      %
      %      d^2 u(x)
      % A E ---------- = -b(x)
      %        dx^2

      dydx = zeros(6,1);
      %         y(1)          = v
      dydx(1) = y(2);       % = theta
      dydx(2) = y(3)/(E*I); % = M/(EI)
      dydx(3) = y(4);       % = V
      dydx(4) = qyloc(x);   % = qyloc
      dydx(5) = y(6)/(A*E); % = u
      dydx(6) = -qxloc(x);  % = faxial
   end

%% ------------------------------------------------------------------------

   function res = condiciones_de_apoyo(YL,YR)
      % condiciones de apoyo (cond. de frontera de la ecuacion diferencial)
      v_  = 1; t_ = 2; M_ = 3; V_ = 4; u_ = 5; fax_ = 6;
      res = [ % YL: apoyo izquierdo (LEFT), YR: apoyo derecho (RIGHT)
              YL(u_) - 0          % uloc(0)     = 0
              YL(v_) - 0          % vloc(0)     = 0
              YL(t_) - 0          % thetaloc(0) = 0
              YR(u_) - 0          % uloc(L)     = 0
              YR(v_) - 0          % vloc(L)     = 0
              YR(t_) - 0 ];       % thetaloc(L) = 0
   end
end
%% end
