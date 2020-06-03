%% Se ejecuta el método anterior, con el proposito de comparar la respuesta
APm1T_TQQL_metodo1
clearvars -except invP_T_metodo1; clc

%        ^ eta (xib2)
%   xib3 |         
%      \ |         
%       \|         
%       (3) (0,1)
%        |\
%        | \
%        5  4
%        |   \
%        |    \
%        |     \
%        |      \
%        |       \            
%       (6)      (5)
%        |         \      
%        |          \       
%        |           \    
%        |            \   
%        6             3    
%        |              \ 
%        |               \
%        |                \
%       (1)--1---(4)---2--(2)------> xi (xib1)
%      (0,0)              (1,0)


%% Se definen las constantes
syms xi eta 
syms gxi1 geta1 gxi2 geta2 gxi3 geta3 gxi4 geta4 gxi5 geta5 gxi6 geta6

XI = 1; ETA = 2;

%% Paso 0: se definen los puntos de colocacion y las direcciones xi_bar_i
% puntos de colocacion
d = sym(sqrt(3))/6;

nod = [ ...
% xi      eta       beta
 0.5-d     0           0    % 1
 0.5+d     0           0    % 2
 0.5+d   0.5-d    3*pi/4    % 3
 0.5-d   0.5+d    3*pi/4    % 4
   0     0.5+d      pi/2    % 5
   0     0.5-d      pi/2 ]; % 6

n = size(nod,1);           % numero de puntos de colocacion
beta = nod(:,3);           % se definen las direcciones
             
%% Paso 1: Se define la matriz A de interpolacion de las deformaciones 
%angulares en el sistema de coordenadas naturales

A = @(xi,eta) [ 1  xi eta   0   0   0 
                 0  0   0   1  xi eta ];   % eq 6.62, 6.116

Ag = [ A(nod(1,XI), nod(1,ETA))
       A(nod(2,XI), nod(2,ETA))     % eq 6.67
       A(nod(3,XI), nod(3,ETA)) 
       A(nod(4,XI), nod(4,ETA))
       A(nod(5,XI), nod(5,ETA)) 
       A(nod(6,XI), nod(6,ETA)) ];             
             
%% Paso 2: se define la matriz T y el vector gpg
TT = cell(1,n);
for i = 1:n
   TT{i} = [ cos(beta(i)) sin(beta(i)) ];
end
T = blkdiag(TT{:})                              % eq 6.65, 6.89

% gamma prima gorrito eq 6.66b
gpg = [ gxi1 geta1 gxi2 geta2 gxi3 geta3 gxi4 geta4 gxi5 geta5 gxi6 geta6 ].';
   
P = T*Ag                                         % eq 6.68

A_invP_T = A(xi,eta)*inv(P)*T

gp = A_invP_T*gpg                                % eq 6.70

%% Se verifica que el metodo 1 y el metodo 2 sean iguales
disp('El error en la próxima resta debe ser despreciable')
double(P\T) - invP_T_metodo1