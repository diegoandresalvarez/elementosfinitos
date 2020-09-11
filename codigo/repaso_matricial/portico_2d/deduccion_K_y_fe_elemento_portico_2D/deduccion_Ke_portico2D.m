% Programa para deducir la matriz de rigidez de un elemento de portico 2D
clear, clc
syms w x L V(x) M(x) t(x) v(x) EI EA

Kflexion = sym(zeros(4));
for i = 1:4
    sol = dsolve(...       
           diff(V,x) == 0,    ... % se definen las ecuaciones diferenciales
           diff(M,x) == V,    ...
           diff(t,x) == M/EI, ... 
           diff(v,x) == t,    ...
           v(0) == (i==1),    ... % con sus respectivas condiciones de 
           t(0) == (i==2),    ... % frontera  
           v(L) == (i==3),    ...           
           t(L) == (i==4));

    Kflexion(:,i) = [ +subs(sol.V, x, 0)       % Yi  se evaluan las 
                      -subs(sol.M, x, 0)       % Mi  reacciones verticales
                      -subs(sol.V, x, L)       % Yj  y los momentos en los
                      +subs(sol.M, x, L) ];    % Mj  apoyos
end
Kaxial = EA/L * [1 -1; -1 1];

K = sym(zeros(6));
K([1 4],[1 4])         = Kaxial;
K([2 3 5 6],[2 3 5 6]) = Kflexion;
K
