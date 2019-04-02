% Programa para encontrar las cargas nodales equivalentes de:

%      ...````.....                                          
%   /|' |  |  |  |  | ..       E, I, A, L constante                                 
%   /|  |  |  |  |  |   ..                                      
%   /|  |  |  |  |  |  | ..                   .......|\
%   /|  |  |  |  |  |  |  |`....         ....` |  |  |\
%   /|  |  |  |  |  |  |  |  |  `......` |  |  |  |  |\
%   /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |\   -q(x) -> variable
%   /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |\
%   /|  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  |\
%   /|###############################################|\
%   /|                                               |\
%
%    |-----------------------L-----------------------|

clear, clc

syms w x L
% syms q;   % si se supone que q es constante y comente la siguiente linea
q = -w*x/L; % carga triangular
syms V(x) M(x) t(x) v(x) EI
sol = dsolve(...       
       diff(V) == q,    ... % se definen las ecuaciones diferenciales
       diff(M) == V,    ...
       diff(t) == M/EI, ... 
       diff(v) == t,    ...
       v(0) == 0,       ... % con sus respectivas condiciones de frontera
       v(L) == 0,       ...
       t(0) == 0,       ...   
       t(L) == 0);
      
disp('q = '); disp(simple(q));
disp('V = '); disp(simple(sol.V));
disp('M = '); disp(simple(sol.M));
disp('t = '); disp(simple(sol.t));
disp('v = '); disp(simple(sol.v));

% Se evaluan las cargas nodales equivalentes
disp('Y1 = '); disp(-subs(sol.V, x, 0));
disp('Y2 = '); disp(+subs(sol.V, x, L));

% Y los momentos nodales equivalentes
disp('M1 = '); disp(+subs(sol.M, x, 0));
disp('M2 = '); disp(-subs(sol.M, x, L));
      