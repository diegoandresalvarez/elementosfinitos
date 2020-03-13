% Programa para deducir la matriz de rigidez de un elemento de viga de EB

clear, clc
syms w x L V(x) M(x) t(x) v(x) EI EA

Kb = sym(zeros(4));
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

    Kb(:,i) = [ +subs(sol.V, x, 0)    % Y1  se evaluan las 
                -subs(sol.M, x, 0)    % M1  reacciones verticales
                -subs(sol.V, x, L)    % Y2  y los momentos en los
                +subs(sol.M, x, L) ]; % M2  apoyos
end
disp('Kb = (EI/L^3)*'); pretty(Kb/(EI/L^3));

