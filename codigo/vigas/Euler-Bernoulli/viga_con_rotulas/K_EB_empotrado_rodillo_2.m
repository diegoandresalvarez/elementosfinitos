% Programa para deducir la matriz de rigidez de un elemento de viga de
% Euler-Bernoulli a partir de la solucion de la ecuacion diferencial
% que tiene una rotula a la derecha

clear, clc
syms x L V(x) M(x) t(x) v(x) EI EA

%% Se calcula la matrix de rigidez para la rotula en el lado derecho
K = sym(zeros(3));
for i = 1:3
    sol = dsolve(...       
           diff(V,x) == 0,    ... % se definen las ecuaciones
           diff(M,x) == V,    ... % diferenciales
           diff(t,x) == M/EI, ... 
           diff(v,x) == t,    ...
           v(0) == (i==1),    ... % con sus respectivas 
           t(0) == (i==2),    ... % condiciones de frontera
           v(L) == (i==3),    ...           
           M(L) == 0);

    K(:,i) = [ +subs(sol.V, x, 0)    % Y1  se evaluan las reacciones
               -subs(sol.M, x, 0)    % M1  verticales y los momentos
               -subs(sol.V, x, L) ]; % Y2  en los apoyos
end
K

%% Se calcula la matriz de rigidez para la rotula en el lado izquierdo
K = sym(zeros(3));
for i = 1:3
    sol = dsolve(...       
           diff(V,x) == 0,    ... % se definen las ecuaciones 
           diff(M,x) == V,    ... % diferenciales
           diff(t,x) == M/EI, ... 
           diff(v,x) == t,    ...
           v(0) == (i==1),    ... % con sus respectivas 
           M(0) == 0,         ... % condiciones de frontera  
           v(L) == (i==2),    ...           
           t(L) == (i==3));

    K(:,i) = [ +subs(sol.V, x, 0)    % Y1  se evaluan las reacciones
               -subs(sol.V, x, L)    % Y2  verticales y los momentos
               +subs(sol.M, x, L) ]; % M2  en los apoyos
end
K
