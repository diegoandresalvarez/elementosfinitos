%% Programa para calcular la matriz de rigidez K de un elemento como el 
%  mostrado
%  /|      E, I, A constante
%  /|________________________________________________
%  /|  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^ q(x) -> variable
%  /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
%  /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
%  /|################################################
%  /|                                               o
%  /|                                              /.\
%  /|                                             /////

%   |-----------------------L-----------------------|

%%
clear, clc

%% Se definen las variables simbolicas y algunas constantes
syms AE EI L u1 v1 t1 u2 v2 t2 
L2 = L^2; 
L3 = L^3;

%% Se define el vector de desplazamientos nodales a y la matriz K para el
% elemento doblemente empotrado

ae = [u1; v1; t1; u2; v2; t2];

Keloc = [  AE/L          0         0  -AE/L         0         0       
              0   12*EI/L3   6*EI/L2      0 -12*EI/L3   6*EI/L2
              0    6*EI/L2   4*EI/L       0  -6*EI/L2   2*EI/L
          -AE/L          0         0   AE/L         0         0
              0  -12*EI/L3  -6*EI/L2      0  12*EI/L3  -6*EI/L2
              0    6*EI/L2   2*EI/L       0  -6*EI/L2   4*EI/L   ];
    
%%  Se multiplica Keloc*ae para obtener las 6 ecuaciones correspondientes    
fe = Keloc*ae;

%% La sexta ecuacion se iguala a cero ya que no hay momentos y se despeja t2
M2 = fe(6);      
t2 = solve(M2 == 0, t2);

%% Se reemplazan los resultados de nuevo en ae y se recalcula Keloc*ae
ae = [u1; v1; t1; u2; v2; t2];
fe = simplify(Keloc*ae);

%% Se elimina el gdl 6 y se recalcula la matriz K:
ee = mat2cell(eye(5), [1 1 1 1 1], 5); % produce {[1 0 0 0 0];[0 1 0 0 0];[0 0 1 0 0];[0 0 0 1 0];[0 0 0 0 1]}
K = subs(fe([1 2 3 4 5]'), {u1, v1, t1, u2, v2}', ee)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Ahora se hace lo mismo pero con:
%        E, I, A constante                           |\
%    ________________________________________________|\
%    ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  ^  |\  q(x) -> variable
%    |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |\
%    |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |\
%    ################################################|\
%    o                                               |\
%   /.\                                              |\
%  /////                                             |\

%    |-----------------------L-----------------------|

clearvars t1 t2 ae

%% Se definen las variables simbolicas y algunas constantes
syms t1 t2 

%% Se define el vector de desplazamientos nodales a para el elemento doblemente empotrado
ae = [u1; v1; t1; u2; v2; t2];

%%  Se multiplica Keloc*ae para obtener las 6 ecuaciones correspondientes    
fe = Keloc*ae;

%% La tercera ecuacion se iguala a cero ya que no hay momentos y se despeja t1
% se despeja t1
M1 = fe(3);      
t1 = solve(M1 == 0, t1);

%% Se reemplazan los resultados de nuevo en ae y se recalcula Keloc*ae
ae = [u1; v1; t1; u2; v2; t2];
fe = simplify(Keloc*ae);

%% Se elimina el gdl 3 y se recalcula la matriz K:
K = subs(fe([1 2 4 5 6]), {u1, v1, u2, v2, t2}', ee)
