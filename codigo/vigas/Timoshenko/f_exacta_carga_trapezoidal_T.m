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
% para una viga de Timoshenko, asumiendo E, I, G, A constante 

clear, clc
syms w x L w1 w2 V(x) M(x) t(x) v(x) EI beta
GAast = (12 * EI)/(L^2 * beta);

q = (w2-w1)*x/L + w1;       % carga trapezoidal

sol = dsolve(...       
       diff(V) == q,    ... % se definen las ecuaciones diferenciales
       diff(M) == V,    ...
       diff(t) == M/EI, ... 
       diff(v) == t - V/GAast,    ...
       v(0) == 0,       ... % con sus respectivas condiciones de frontera
       v(L) == 0,       ... % para una viga doblemente empotrada
       t(0) == 0,       ...   
       t(L) == 0);
      
disp('q = '); disp(simplify(q));
disp('V = '); disp(simplify(sol.V));
disp('M = '); disp(simplify(sol.M));
disp('t = '); disp(simplify(sol.t));
disp('v = '); disp(simplify(sol.v));

% Se evaluan las cargas nodales equivalentes
Y1 = simplify(-subs(sol.V, x, 0));
Y2 = simplify(+subs(sol.V, x, L));
disp('Y1 = '); pretty(Y1);
disp('Y2 = '); pretty(Y2);

% Y los momentos nodales equivalentes
M1 = simplify(+subs(sol.M, x, 0));
M2 = simplify(-subs(sol.M, x, L));
disp('M1 = '); pretty(M1);
disp('M2 = '); pretty(M2);

%% Se imprime el vector de fuerzas nodales equivalentes
f = [Y1; M1; Y2; M2];
coef = 1/(beta + 1);
disp('f_TE = 1/(beta + 1) *')
pretty(simplify(expand(f/coef)));

%% Se compara dicho vector con la teoria de Euler-Bernoulli
disp('f_EB = ')
pretty(simplify(limit(f, beta, 0)))