APm1T_QLLL_metodo1

clearvars -except A_invP_T_metodo1 gp_metodo1; 
clc

%          ^ eta
%          |
%          |
%          3
%   (4)----+----(3)
%    |           |  
%    |           |
%  4 x           x 2  ------> xi
%    |           |
%    |           |
%   (1)----+----(2)
%          1

%% Se definen las constantes
syms xi eta gxi1 geta1 gxi2 geta2 gxi3 geta3 gxi4 geta4

XI = 1; ETA = 2;

%% Paso 0: se definen los puntos de colocacion y las direcciones xi_bar_i
% puntos de colocacion
nod = [...
% xi  eta    beta_i
   0   -1    0         % 1 +
   1    0    pi/2      % 2 x
   0    1    pi        % 3 + 
  -1    0    pi/2 ];   % 4 x
          
% se definen las direcciones          
beta = sym(nod(:,3));

%% Paso 1: Se define la matriz A de interpolacion de las deformaciones 
%angulares en el sistema de coordenadas naturales

A = @(xi,eta) [ 1  eta   0    0
                0    0   1   xi ];         % eq 6.62, 6.68, 6.86

Ag = [ A(nod(1,XI), nod(1,ETA))
       A(nod(2,XI), nod(2,ETA))            % eq 6.67
       A(nod(3,XI), nod(3,ETA)) 
       A(nod(4,XI), nod(4,ETA)) ];             
             
%% Paso 2: se define la matriz T y el vector gpg
ngamma = length(beta);                     % numero de puntos de colocacion
TT = cell(1,ngamma);
for i = 1:ngamma
   TT{i} = [cos(beta(i)) sin(beta(i)) ];
end
T = blkdiag(TT{:})                         % eq 6.65, 6.89

% gamma prima gorrito eq 6.66b
gpg = [gxi1; geta1; gxi2; geta2; gxi3; geta3; gxi4; geta4];
   
P = T*Ag                                   % eq 6.68
A_invP_T = A(xi,eta)*inv(P)*T
gp = A_invP_T*gpg                          % eq 6.70

% se verifica que el metodo 1 y 2 dan los mismos resultados
disp(simplify(gp - gp_metodo1))
disp(simplify(A_invP_T - A_invP_T_metodo1))
