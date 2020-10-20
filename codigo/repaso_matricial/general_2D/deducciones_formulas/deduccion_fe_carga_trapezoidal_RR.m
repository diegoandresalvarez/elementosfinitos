% Programa para encontrar las cargas nodales equivalentes de:

%  -q1 
%   |  -  _
%   |  |  |  -  _ 
%   |  |  |  |  |  -  _                   E, I, A constante
%   |  |  |  |  |  |  |  -  _
%   |  |  |  |  |  |  |  |  |  -  _
%   |  |  |  |  |  |  |  |  |  |  | -q2  q(x) -> carga vertical trapezoidal
%   |  |  |  |  |  |  |  |  |  |  |  |
%   V  V  V  V  V  V  V  V  V  V  V  V
%   b1                              b2   b(x) -> carga axial trapezoidal
%    ----> ----> ---> ---> --> --> -> 
%   o################################o
%  / \                              / \
% -----                            -----
% /////                            /////
%   |                                |
%   x1                               x2
%   |---------------L----------------|

clear, clc

syms x L q1 q2 b1 b2 m b
x1 = 0;
x2 = L;
sol = solve(q1 == m*x1+b, q2 == m*x2+b, m, b);
q = simplify(sol.m*x + sol.b); % carga vertical trapezoidal

sol = solve(b1 == m*x1+b, b2 == m*x2+b, m, b);
b = simplify(sol.m*x + sol.b); % carga axial trapezoidal

syms fax(x) V(x) M(x) t(x) u(x) v(x) E A I u1 u2 v1 v2
EI = E*I;
EA = E*A;
sol = dsolve(...       
       diff(V) == q,      ... % se definen las ecuaciones diferenciales
       diff(M) == V,      ...
       diff(t) == M/EI,   ... 
       diff(v) == t,      ...
       diff(u) == fax/EA, ...
       diff(fax) == -b,   ...
       u(0) == 0,         ...
       u(L) == 0,         ...
       v(0) == 0,         ... % con sus respectivas condiciones de frontera
       v(L) == 0,         ...
       M(0) == 0,         ...   
       M(L) == 0);

%{   
disp('q = ');   disp(q);
disp('fax = '); disp(simplify(sol.fax));
disp('V = ');   disp(simplify(sol.V));
disp('M = ');   disp(simplify(sol.M));
disp('t = ');   disp(simplify(sol.t));
disp('u = ');   disp(simplify(sol.u));
disp('v = ');   disp(simplify(sol.v));
%}

% Se evaluan las cargas nodales equivalentes
disp('X1 = '); disp(simplify(+subs(sol.fax, x, 0)));
disp('X2 = '); disp(simplify(-subs(sol.fax, x, L)));

disp('Y1 = '); disp(simplify(-subs(sol.V, x, 0)));
disp('Y2 = '); disp(simplify(+subs(sol.V, x, L)));

% Y los momentos nodales equivalentes
disp('M1 = '); disp(simplify(+subs(sol.M, x, 0)));
disp('M2 = '); disp(simplify(-subs(sol.M, x, L)));
      
disp('-------------------------------------------------------------------')

sol = dsolve(...       
       diff(V) == q,      ... % se definen las ecuaciones diferenciales
       diff(M) == V,      ...
       diff(t) == M/EI,   ... 
       diff(v) == t,      ...
       diff(u) == fax/EA, ...
       diff(fax) == -b,   ...
       u(0) == u1,         ...
       u(L) == u2,         ...
       v(0) == v1,         ... % con sus respectivas condiciones de frontera
       v(L) == v2,         ...
       M(0) == 0,         ...   
       M(L) == 0);
      
disp('q = ');   disp(q);
disp('fax = '); disp(simplify(sol.fax));
disp('V = ');   disp(simplify(sol.V));
disp('M = ');   disp(simplify(sol.M));
disp('t = ');   disp(simplify(sol.t));
disp('u = ');   disp(simplify(sol.u));
disp('v = ');   disp(simplify(sol.v));

%{
% Se evaluan las cargas nodales equivalentes
disp('X1 = '); disp(simplify(+subs(sol.fax, x, 0)));
disp('X2 = '); disp(simplify(-subs(sol.fax, x, L)));

disp('Y1 = '); disp(simplify(-subs(sol.V, x, 0)));
disp('Y2 = '); disp(simplify(+subs(sol.V, x, L)));

% Y los momentos nodales equivalentes
disp('M1 = '); disp(simplify(+subs(sol.M, x, 0)));
disp('M2 = '); disp(simplify(-subs(sol.M, x, L)));
%}      