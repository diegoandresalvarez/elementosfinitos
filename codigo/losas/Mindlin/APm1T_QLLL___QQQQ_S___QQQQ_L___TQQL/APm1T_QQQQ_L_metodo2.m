%% Se ejecuta el método anterior, con el proposito de comparar la respuesta
APm1T_QQQQ_L_metodo1
clearvars -except A_invP_T_metodo1 gp_metodo1; clc

%             ^ eta
%             |
%             |
%             |
%   (7)---2--(6)--1---(5)
%    |                 |  
%    |                 |
%   11        9        7  
%    |                 |
%   (8)   4  (9)  3   (4) ------> xi
%    |                 |  
%   12       10        8
%    |                 |
%    |                 |
%   (1)---6--(2)--5---(3)

%% Se definen las constantes
syms xi eta gxi1 geta1 gxi2 geta2 gxi3 geta3 gxi4 geta4 gxi5 geta5 gxi6 geta6 
syms        gxi7 geta7 gxi8 geta8 gxi9 geta9 gxi10 geta10 gxi11 geta11 gxi12 geta12

XI = 1; ETA = 2;

%% Paso 0: se definen los puntos de colocacion y las direcciones xi_bar_i
% puntos de colocacion
%syms a
a = 1/sym(sqrt(3));
           
            nod = [ ...
% xi eta
   a   1     pi      %  1 +
  -a   1     pi      %  2 +
   a   0      0      %  3
  -a   0      0      %  4
   a  -1      0      %  5 +
  -a  -1      0      %  6 +
   1   a   pi/2      %  7 x
   1  -a   pi/2      %  8 x        
   0   a   pi/2      %  9 
   0  -a   pi/2      % 10        
  -1   a   pi/2      % 11 x
  -1  -a   pi/2 ];   % 12 x
          
% se definen las direcciones          
beta = sym(nod(:,3));    % angulos
n = length(beta);        % numero puntos de colocacion

%% Paso 1: Se define la matriz A de interpolacion de las deformaciones 
%angulares en el sistema de coordenadas naturales

A = @(xi,eta) [ 1  xi  eta  xi*eta  eta^2  xi*eta^2   0   0    0       0     0        0
                 0  0    0       0      0         0   1  xi  eta  xi*eta  xi^2 eta*xi^2 ];   % eq 6.62, 6.99

Ag = [ A(nod( 1,XI), nod( 1,ETA))
       A(nod( 2,XI), nod( 2,ETA))     % eq 6.67
       A(nod( 3,XI), nod( 3,ETA)) 
       A(nod( 4,XI), nod( 4,ETA))
       A(nod( 5,XI), nod( 5,ETA)) 
       A(nod( 6,XI), nod( 6,ETA)) 
       A(nod( 7,XI), nod( 7,ETA)) 
       A(nod( 8,XI), nod( 8,ETA))
       A(nod( 9,XI), nod( 9,ETA))
       A(nod(10,XI), nod(10,ETA))
       A(nod(11,XI), nod(11,ETA))       
       A(nod(12,XI), nod(12,ETA)) ];             
             
%% Paso 2: se define la matriz T y el vector gpg
TT = cell(1,n);
for i = 1:n
   TT{i} = [ cos(beta(i)) sin(beta(i)) ];
end
T = blkdiag(TT{:})                              % eq 6.65, 6.89

% gamma prima gorrito eq 6.66b
gpg = [ gxi1 geta1 gxi2 geta2 gxi3 geta3 gxi4 geta4 gxi5 geta5 ...
        gxi6 geta6 gxi7 geta7 gxi8 geta8 gxi9 geta9 gxi10 geta10 ...
        gxi11 geta11 gxi12 geta12 ].';     
   
P = T*Ag                                         % eq 6.68

A_invP_T = factor(A(xi,eta)*inv(P)*T)

gp = A_invP_T*gpg                                % eq 6.70

disp(simple(gp-gp_metodo1))

disp(simple(A_invP_T - A_invP_T_metodo1))

disp(simple((A_invP_T- A_invP_T_metodo1)*gpg))
