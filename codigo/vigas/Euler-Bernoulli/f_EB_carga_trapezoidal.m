% Programa para encontrar las cargas nodales equivalentes asociados a una carga
% trapezoidal:
%
%   /|                                      ^  q(x)
%   /|                             ^  |  |  |                   
%   /|                    ^  |  |  |  |  |  |  |/
%   /|           ^  |  |  |  |  |  |  |  |  |  |/
%   /|  ^  |  |  |  |  |  |  |  |  |  |  |  |  |/
%   /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |/  
%   /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |/
%   /|  |  |  |  |  |  |  |  |  |  |  |  |  |  |/
%   /|#########################################|/
%   /|                                         |/
%
%    |--------------------L--------------------|
% 
% para una viga de Euler-Bernoulli, asumiendo E, I, A constante 

% FECHA         QUIEN  QUE 
% Ago 12, 2022  DAAM   El código es igual que el de PYTHON
%
% DAAM >>> Diego Andrés Alvarez Marín daalvarez@unal.edu.co

clear, clc

syms syms w(x) x L w1 w2 V(x) M(x) t(x) EI

q = (w2-w1)*x/L + w1;       % carga trapezoidal

sol = dsolve(...       
       diff(V,x) == q,    ... % se definen las ecuaciones diferenciales
       diff(M,x) == V,    ...
       diff(t,x) == M/EI, ... 
       diff(w,x) == t,    ...
       w(0) == 0,         ... % con sus respectivas condiciones de frontera
       w(L) == 0,          ...
       t(0) == 0,         ...   
       t(L) == 0);

% se imprimen los resultados
disp('q = '); pretty(simplify(q));
disp('V = '); pretty(simplify(sol.V));
disp('M = '); pretty(simplify(sol.M));
disp('t = '); pretty(simplify(sol.t));
disp('w = '); pretty(simplify(sol.w));

% Se evaluan las cargas nodales equivalentes
disp('Y1 = '); pretty(simplify(-subs(sol.V, x, 0)));
disp('Y2 = '); pretty(simplify(+subs(sol.V, x, L)));

% Y los momentos nodales equivalentes
disp('M1 = '); pretty(simplify(+subs(sol.M, x, 0)));
disp('M2 = '); pretty(simplify(-subs(sol.M, x, L)));
      
