% Programa para deducir la matriz de rigidez de un elemento de viga de
% Euler-Bernoulli a partir de la solucion de la ecuacion diferencial
clear, clc
syms x L V(x) M(x) t(x) v(x) EI EA

%% Se calcula la matrix de rigidez para la rotula en el lado derecho
K = sym(zeros(3));
for i = 1:3
    sol = dsolve(...       
           diff(V,x) == 0,    ... % se definen las ecuaciones diferenciales
           diff(M,x) == V,    ...
           diff(t,x) == M/EI, ... 
           diff(v,x) == t,    ...
           v(0) == (i==1),    ... % con sus respectivas condiciones de 
           t(0) == (i==2),    ... % frontera  
           v(L) == (i==3),    ...           
           M(L) == 0);

    K(:,i) = [ +subs(sol.V, x, 0)    % Y1  se evaluan las 
               -subs(sol.M, x, 0)    % M1  reacciones verticales
               -subs(sol.V, x, L) ]; % Y2  y los momentos en los apoyos
end
K

%% Se calcula la matrix de rigidez para la rotula en el lado izquierdo
K = sym(zeros(3));
for i = 1:3
    sol = dsolve(...       
           diff(V,x) == 0,    ... % se definen las ecuaciones diferenciales
           diff(M,x) == V,    ...
           diff(t,x) == M/EI, ... 
           diff(v,x) == t,    ...
           v(0) == (i==1),    ... % con sus respectivas condiciones de 
           M(0) == 0,         ... % frontera  
           v(L) == (i==2),    ...           
           t(L) == (i==3));

    K(:,i) = [ +subs(sol.V, x, 0)    % Y1  se evaluan las 
               -subs(sol.V, x, L)    % Y2  reacciones verticales
               +subs(sol.M, x, L) ]; % M2  y los momentos en los apoyos
end
K
