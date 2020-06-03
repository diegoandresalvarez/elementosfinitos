%% Se ejecuta el método anterior, con el proposito de comparar la respuesta
APm1T_QQQQ_S_metodo1
clearvars -except A_invP_T_metodo1 gp_metodo1; clc

%             ^ eta
%             |
%             |
%             |
%   (7)---2--(6)--1---(5)
%    |                 |  
%    |                 |
%    8                 6  
%    |                 |
%   (8)       3       (4) ------> xi
%    |                 |  
%    9                 7
%    |                 |
%    |                 |
%   (1)---5--(2)--4---(3)

%% Se definen las constantes
syms xi eta gxi1 geta1 gxi2 geta2 gxi3 geta3 gxi4 geta4
syms        gxi5 geta5 gxi6 geta6 gxi7 geta7 gxi8 geta8 gxi9 geta9

XI = 1; ETA = 2;

%% Paso 0: se definen los puntos de colocacion y las direcciones xi_bar_i
% puntos de colocacion
%syms a
a = 1/sym(sqrt(3));
           
nod = [ ...
% xi eta
   a   1    % 1
  -a   1    % 2
   0   0    % 3
   a  -1    % 4
  -a  -1    % 5
   1   a    % 6
   1  -a    % 7              
  -1   a    % 8
  -1  -a ]; % 9
          
% se definen las direcciones          
theta = sym([     pi      % 1
                  pi      % 2
                  pi      % 3h
                pi/2      % 3v
                   0      % 4
                   0      % 5
                pi/2      % 6 
                pi/2      % 7                 
                pi/2      % 8
                pi/2 ]);  % 9 
             
%% Paso 1: Se define la matriz A de interpolacion de las deformaciones 
%angulares en el sistema de coordenadas naturales

A = @(xi,eta) [ 1  xi  eta  xi*eta  eta^2  0   0    0       0     0 
                 0  0    0       0      0  1  xi  eta  xi*eta  xi^2 ];   % eq 6.62, 6.94

Ag = [ A(nod(1,XI), nod(1,ETA))
       A(nod(2,XI), nod(2,ETA))     % eq 6.67
       A(nod(3,XI), nod(3,ETA)) 
       A(nod(3,XI), nod(3,ETA)) 
       A(nod(4,XI), nod(4,ETA))
       A(nod(5,XI), nod(5,ETA)) 
       A(nod(6,XI), nod(6,ETA)) 
       A(nod(7,XI), nod(7,ETA)) 
       A(nod(8,XI), nod(8,ETA))
       A(nod(9,XI), nod(9,ETA)) ];             
             
%% Paso 2: se define la matriz T y el vector gpg
ntheta = length(theta);
TT = cell(1,ntheta);
for i = 1:ntheta
   TT{i} = [ cos(theta(i)) sin(theta(i)) ];
end
T = blkdiag(TT{:})                              % eq 6.65, 6.89

% gamma prima gorrito eq 6.66b
gpg = [ gxi1 geta1 gxi2 geta2 gxi3 geta3 gxi3 geta3 gxi4 geta4 ...
        gxi5 geta5 gxi6 geta6 gxi7 geta7 gxi8 geta8 gxi9 geta9 ].';
   
P = T*Ag                                         % eq 6.68

A_invP_T = A(xi,eta)*inv(P)*T

gp = A_invP_T*gpg                                % eq 6.70


%% Se compacta la matriz A_invP_T (se remueve la repeticion de "gxi3 geta3")
A_invP_T(:,5) = A_invP_T(:,5) + A_invP_T(:,7);
A_invP_T(:,6) = A_invP_T(:,6) + A_invP_T(:,8);
A_invP_T(:,[7 8]) = [];

%% Se verifica que el metodo 1 y el metodo 2 sean iguales
gpg = [ gxi1 geta1 gxi2 geta2 gxi3 geta3            gxi4 geta4 ...
        gxi5 geta5 gxi6 geta6 gxi7 geta7 gxi8 geta8 gxi9 geta9 ].';

disp(simple(gp-gp_metodo1))

disp(simple(A_invP_T - A_invP_T_metodo1))

disp(simple((A_invP_T- A_invP_T_metodo1)*gpg))
