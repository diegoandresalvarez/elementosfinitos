% Programa para deducir la matriz de rigidez de un elemento de viga 2D de
% Timoshenko a partir de la solucion de la ecuacion diferencial
clear, clc
syms q x L V(x) M(x) t(x) v(x) EI EA GAast

%% Se calcula la matrix de rigidez
K_T = sym(zeros(4));
Nv = sym(zeros(1,4)); % func forma para los desplazamientos
Nt = sym(zeros(1,4)); % func forma para los giros de la seccion transversal
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

    K_T(:,i) = [ +subs(sol.V, x, 0)    % Y1  se evaluan las 
                 -subs(sol.M, x, 0)    % M1  reacciones verticales
                 -subs(sol.V, x, L)    % Y2  y los momentos en los
                 +subs(sol.M, x, L) ]; % M2  apoyos

    Nv(i) = sol.v; % v = Nv * [v1; t1 v2 t2]
    Nt(i) = sol.t; % t = Nt * [v1; t1 v2 t2]             
end

%% Se imprime la solucion
beta = (12 * EI)/(L^2 * GAast);
tmp = EI/((1 + beta)*L^3);
disp('K_T = (EI/((1 + beta)*L^3)) * ');
pretty(simplify(K_T/tmp))
% Nota: se puede demostrar que:
% K22 = K44 = simplify((4+beta)*L^2)
% K24 = K42 = simplify((2-beta)*L^2)

%{
% Comprobacion de que a partir de las funciones de forma se puede obtener
% la matriz de rigidez K_T
Bb = diff(Nt, x);
Bs = diff(Nv, x) - Nt;
K_T2 = int(Bb.'*EI*Bb,x,0,L) + int(Bs.'*GAast*Bs,x,0,L);
simplify(K_T - K_T2)
%}

%% Se calcula la matriz de rigidez de Euler-Bernoulli
% Observe que cuando GAast -> Inf, la matriz de rigidez K se vuelve la 
% misma matriz de rigidez K de la teoria de Euler-Bernoulli:
K_EB = limit(K_T, GAast, inf);
disp('K_EB = (EI/L^3) * ');
pretty(simplify(K_EB/(EI/L^3)))

%% Se calculan las fuerzas nodales equivalentes:
syms w
% q = w;
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

f_T = [ -subs(sol.V, x, 0)   % Yi  se evaluan las reacciones verticales y los
        +subs(sol.M, x, 0)   % Mi  momentos en los apoyos y se les multiplica
        +subs(sol.V, x, L)   % Yj  por -1 para estimar la fuerza nodal 
        -subs(sol.M, x, L) ];% Mj  equivalente

%% Y en el limitem cuando beta -> inf, obtenemos el vector f para EB
f_EB = limit(f_T, GAast, inf);
disp('f_EB =');
pretty(simplify(f_EB))
    
%%    
clear beta
syms beta
f_T = simplify(subs(f_T, GAast, (12 * EI)/(L^2 * beta)));
disp('f_T = ');
pretty(f_T)

%% bye, bye!