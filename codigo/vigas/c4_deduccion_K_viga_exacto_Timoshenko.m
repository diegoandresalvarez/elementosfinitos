% Programa para deducir la matriz de rigidez de un elemento de viga 2D de
% Timoshenko a partir de la solucion de la ecuacion diferenfial
clear, clc
syms q w x L V(x) M(x) t(x) v(x) EI EA GAast

%% Se calcula la matrix de rigidez
K = sym(zeros(4));
for i = 1:4
    sol = dsolve(...       
           diff(V,x) == 0,    ... % se definen las ecuaciones diferenciales
           diff(M,x) == V,    ...
           diff(t,x) == M/EI, ... 
           diff(v,x) == t - V/GAast,    ...
           v(0) == (i==1),    ... % con sus respectivas condiciones de 
           t(0) == (i==2),    ... % frontera  
           v(L) == (i==3),    ...           
           t(L) == (i==4));

    K(i,:) = [ +subs(sol.V, x, 0), ...  % Yi  se evaluan las 
               -subs(sol.M, x, 0), ...  % Mi  reacciones verticales
               -subs(sol.V, x, L), ...  % Yj  y los momentos en los
               +subs(sol.M, x, L) ];    % Mj  apoyos
end
beta = (12 * EI)/(L^2 * GAast);
coef1 = EI/((1 + beta)*L^3);
disp('K = (EI/((1 + beta)*L^3)) * ');
pretty(simplify(K/coef1))
% Nota: se puede demostrar que:
% K22 = K44 = simplify((4+beta)*L^2)
% K24 = K42 = simplify((2-beta)*L^2)

%% Se calculan las fuerzas nodales equivalentes:
q = -w*x/L;
sol = dsolve(...       
       diff(V,x) == q,    ... % se definen las ecuaciones diferenciales
       diff(M,x) == V,    ...
       diff(t,x) == M/EI, ... 
       diff(v,x) == t - V/GAast, ...
       v(0) == 0,         ... % condiciones de frontera para una viga 
       t(0) == 0,         ... % doblemente empotrada
       v(L) == 0,         ...           
       t(L) == 0); 

f = [ -subs(sol.V, x, 0)   % Yi  se evaluan las reacciones verticales y los
      +subs(sol.M, x, 0)   % Mi  momentos en los apoyos y se les multiplica
      +subs(sol.V, x, L)   % Yj  por -1 para estimar la fuerza nodal 
      -subs(sol.M, x, L) ];% Mj  equivalente
coef2 = GAast*L^2  + 12*EI;
disp('f = 1/(GAast*L^2  + 12*EI) * ');
pretty(simplify(f*coef2))
