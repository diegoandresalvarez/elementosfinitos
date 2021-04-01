% This is a demostration for equation 24 (in paper Part II) in matrix form

clear
clc

%          ^ eta
%          |
%          |
%          7
%   (4)----+----(3)
%    |           |  
%    |           |
%  8 x           x 6  ------> xi
%    |           |
%    |           |
%   (1)----+----(2)
%          5

%% Se definen las constantes
syms xi eta gxi5 geta5 gxi6 geta6 gxi7 geta7 gxi8 geta8

XI = 1; ETA = 2;

%% Paso 0: se definen los puntos de colocacion y las direcciones xi_bar_i
% puntos de colocacion
nod = [...
% xi  eta    beta_i
   0   -1       0      % 5 +
   1    0    pi/2      % 6 x
   0    1       0      % 7 + 
  -1    0    pi/2 ];   % 8 x
          
% se definen las direcciones
beta = sym(nod(:,3));

%% Paso 1: Se define la matriz A de interpolacion de las deformaciones 
%angulares en el sistema de coordenadas naturales

A = @(xi,eta) [ 1  eta   0    0
                0    0   1   xi ];         % eq 6.62, 6.68, 6.86

Ag = [ A(nod(5-4,XI), nod(5-4,ETA))
       A(nod(6-4,XI), nod(6-4,ETA))            % eq 6.67
       A(nod(7-4,XI), nod(7-4,ETA)) 
       A(nod(8-4,XI), nod(8-4,ETA)) ];             
             
%% Paso 2: se define la matriz T y el vector gpg
ngamma = length(beta);                     % numero de puntos de colocacion
TT = cell(1,ngamma);
for i = 1:ngamma
   TT{i} = [cos(beta(i)) sin(beta(i)) ];
end
T = blkdiag(TT{:});                        % eq 6.65, 6.89

% gamma prima gorrito eq 6.66b
gpg = [gxi5; geta5; gxi6; geta6; gxi7; geta7; gxi8; geta8];
   
P = T*Ag                                   % eq 6.68
A_invP_T = A(xi,eta)*inv(P)*T

% And this is equation 24
gp = A_invP_T*gpg                          % eq 6.70 // eq. 24
