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

%% Se definen las variables simbólicas y algunas constantes
syms EI L v1 t1 v2 t2 
L2 = L^2; 
L3 = L^3;

%% Se define el vector de desplazamientos nodales a y la matriz K para el
% elemento doblemente empotrado

ae = [v1; t1; v2; t2];

Keloc = [  12*EI/L3   6*EI/L2  -12*EI/L3   6*EI/L2
            6*EI/L2   4*EI/L    -6*EI/L2   2*EI/L
          -12*EI/L3  -6*EI/L2   12*EI/L3  -6*EI/L2
            6*EI/L2   2*EI/L    -6*EI/L2   4*EI/L   ];
    
%%  Se multiplica Keloc*ae para obtener las 4 ecuaciones correspondientes    
fe = Keloc*ae;

%% La cuarta ecuacion se iguala a cero ya que no hay momentos y se despeja t2
M2 = fe(4);      
t2 = solve(M2 == 0, t2);

%% Se reemplazan los resultados de nuevo en ae y se recalcula Keloc*ae
ae = [v1; t1; v2; t2];
fe = simplify(Keloc*ae);

%% Se elimina el gdl 4 (el de M2) y se recalcula la matriz K:
ee = mat2cell(eye(3), [1 1 1], 3); % produce {[1 0 0];[0 1 0];[0 0 1]}
K = subs(fe([1 2 3]'), {v1, t1, v2}', ee)

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

%% Se definen las variables simbólicas y algunas constantes
syms t1 t2 

%% Se define el vector de desplazamientos nodales a para el elemento doblemente empotrado
ae = [v1; t1; v2; t2];

%%  Se multiplica Keloc*ae para obtener las 4 ecuaciones correspondientes    
fe = Keloc*ae;

%% La segunda ecuacion se iguala a cero ya que no hay momentos y se despeja t1
% se despeja t1
M1 = fe(2);      
t1 = solve(M1 == 0, t1);

%% Se reemplazan los resultados de nuevo en ae y se recalcula Keloc*ae
ae = [v1; t1; v2; t2];
fe = simplify(Keloc*ae);

%% Se elimina el gdl 2 (el de M1) y se recalcula la matriz K:
K = subs(fe([1 3 4]), {v1, v2, t2}', ee)
