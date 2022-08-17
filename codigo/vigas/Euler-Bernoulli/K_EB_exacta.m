% Programa para deducir la matriz de rigidez de un elemento de viga de
% Euler-Bernoulli a partir de la solucion de la ecuacion diferencial

% FECHA         QUIEN  QUE 
% Ago 11, 2022  DAAM   El código es igual que el de PYTHON
%
% DAAM >>> Diego Andrés Alvarez Marín daalvarez@unal.edu.co

clear, clc
syms x L V(x) M(x) t(x) w(x) EI

%% Se calcula la matrix de rigidez
K = sym(zeros(4));
for i = 1:4
    sol = dsolve(...       
           diff(V,x) == 0,    ... % se definen las ecuaciones diferenciales
           diff(M,x) == V,    ...
           diff(t,x) == M/EI, ... 
           diff(w,x) == t,    ...
           w(0) == (i==1),    ... % con sus respectivas condiciones de 
           t(0) == (i==2),    ... % frontera  
           w(L) == (i==3),    ...           
           t(L) == (i==4));

    K(:,i) = [ +subs(sol.V, x, 0)    % Y1  se evaluan las 
               -subs(sol.M, x, 0)    % M1  reacciones verticales
               -subs(sol.V, x, L)    % Y2  y los momentos en los
               +subs(sol.M, x, L) ]; % M2  apoyos
end

%% Se imprime la solucion
disp('K_EB = (EI/L^3) *'); pretty(K/(EI/L^3));

%% Bye, bye!