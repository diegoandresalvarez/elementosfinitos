% Programa para deducir la matriz de rigidez de un elemento de viga 2D de
% Timoshenko a partir de la solucion de la ecuacion diferencial
clear, clc, close all
syms q x L V(x) M(x) t(x) v(x) EI EA beta
GAast = (12 * EI)/(L^2 * beta);

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

    Nv(i) = sol.v; % v = Nv * [v1; t1; v2; t2]
    Nt(i) = sol.t; % t = Nt * [v1; t1; v2; t2]             
end

%% Se imprime la solucion
tmp = EI/(L^3 * (beta + 1));
disp('K_T = EI/(L^3 * (beta + 1)) * ');
pretty(simplify(K_T/tmp))

%%
%{
% Comprobacion de que a partir de las funciones de forma se puede obtener
% la matriz de rigidez K_T
Bb = diff(Nt, x);
Bs = diff(Nv, x) - Nt;
K_T2 = int(Bb.'*EI*Bb,x,0,L) + int(Bs.'*GAast*Bs,x,0,L);
simplify(K_T - K_T2)
%}

%% Se calcula la matriz de rigidez de Euler-Bernoulli
% Observe que cuando GAast -> Inf (o alternativamente beta -> 0), la matriz 
% de rigidez K se vuelve la misma matriz de rigidez K de la teoria de 
% Euler-Bernoulli:
K_EB = limit(K_T, beta, 0);
disp('K_EB = (EI/L^3) * ');
pretty(simplify(K_EB/(EI/L^3)))

%% se calculan las funciones de forma en funcion de xi
syms xi
disp('Se imprimen las funciones de forma')
Nv = simplify(expand(subs(Nv, x, L*(1+xi)/2))).';
Nt = simplify(subs(Nt, x, L*(1+xi)/2)).';
coef = (L^3 * (beta + 1));
disp('Nw^T = 1/(8 * L^3 * (beta + 1)) * '); pretty(collect(8*Nv*coef, xi));
disp('Nt^T = 1/(4 * L^3 * (beta + 1)) * '); pretty(collect(4*Nt*coef, xi));

% se grafican
Nv2 = expand(subs(Nv, {L, beta}, {1, 1}));
figure; fplot(Nv2, [-1,1], 'LineWidth', 2); 
title('Funciones de forma Nw(x) (para L = 1, \beta = 1)')
legend('Nw1(x)','Nw2(x)','Nw3(x)','Nw4(x)')

Nt2 = expand(subs(Nt, {L, beta}, {1, 1}));
figure; fplot(Nt2, [-1,1], 'LineWidth', 2); 
title('Funciones de forma Nt(x) (para L = 1, \beta = 1)')
legend('Nt1(x)','Nt2(x)','Nt3(x)','Nt4(x)')

%% observe que cuando beta -> 0, Nv contiene de las funciones de forma de 
% la viga de Euler-Bernoulli
N_EB = expand(subs(Nv, beta, 0));
disp('N_EB = '); pretty(N_EB)

% Y observe la relacion entre las derivadas
dNEB_dx1 = expand(diff(N_EB*2/L, xi));
dNEB_dx2 = expand(subs(Nt, beta, 0));
dNEB_dx1 - dNEB_dx2

%% Calculo de la matriz de masa consistente "exacta"
syms rho A
NN = Nv.';
M = simplify(int(rho*A*(NN.')*NN*L/2,xi,-1,1))

%% Bye, bye!
return;