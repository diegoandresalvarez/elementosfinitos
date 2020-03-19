% Programa para encontrar las cargas nodales equivalentes asociados a una carga
% trapezoidal:
%
%   /|                                      ^  q(x)
%   /|                             ^  |  |  |                   
%   /|                    ^  |  |  |  |  |  |  |\
%   /|           ^  |  |  |  |  |  |  |  |  |  |\
%   /|  ^  |  |  |  |  |  |  |  |  |  |  |  |  |\
%   /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |\  
%   /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |\
%   /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |\
%   /|#########################################|\
%   /|                                         |\
%
%    |--------------------L--------------------|
% 
% para una viga de Euler-Bernoulli, asumiendo E, I, A constante 

clear, clc

syms syms w x L w1 w2 V(x) M(x) t(x) v(x) EI

q = (w2-w1)*x/L + w1;       % carga trapezoidal

sol = dsolve(...       
       diff(V) == q,    ... % se definen las ecuaciones diferenciales
       diff(M) == V,    ...
       diff(t) == M/EI, ... 
       diff(v) == t,    ...
       v(0) == 0,       ... % con sus respectivas condiciones de frontera
       v(L) == 0,       ...
       t(0) == 0,       ...   
       t(L) == 0);
      
disp('q = '); disp(simplify(q));
disp('V = '); disp(simplify(sol.V));
disp('M = '); disp(simplify(sol.M));
disp('t = '); disp(simplify(sol.t));
disp('v = '); disp(simplify(sol.v));

% Se evaluan las cargas nodales equivalentes
disp('Y1 = '); disp(simplify(-subs(sol.V, x, 0)));
disp('Y2 = '); disp(simplify(+subs(sol.V, x, L)));

% Y los momentos nodales equivalentes
disp('M1 = '); disp(simplify(+subs(sol.M, x, 0)));
disp('M2 = '); disp(simplify(-subs(sol.M, x, L)));
      
