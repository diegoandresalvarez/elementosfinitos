% Programa para deducir la matriz de rigidez de un elemento de viga de
% Euler-Bernoulli+Winkler a partir de la solución de la ecuación diferencial

% FECHA         QUIEN  QUE 
% Ago 14, 2022  DAAM   Código inicial
%
% DAAM >>> Diego Andrés Alvarez Marín daalvarez@unal.edu.co

clear, clc

%% Método 1:
syms x V(x) M(x) t(x) w(x) EI L k

% Se calcula la matrix de rigidez
K = sym(zeros(4));
for i = 1:4
    sol = dsolve(...       
          diff(V,x) + k*w == 0,    ... % se definen las ecuaciones diferenciales
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

K1 = K;

%% Método 2:
lambda = L*(k/(4*EI))^(1/4);
den = (sinh(lambda)^2 - sin(lambda)^2);
k11 = (4*EI/L^3)*lambda^3*(sinh(lambda)*cosh(lambda) +  sin(lambda)* cos(lambda))/den;
k12 = (2*EI/L^2)*lambda^2*(sinh(lambda)^2            +  sin(lambda)^2           )/den;
k13 = (4*EI/L^3)*lambda^3*(sinh(lambda)* cos(lambda) +  sin(lambda)*cosh(lambda))/den;
k14 = (4*EI/L^2)*lambda^2*(sinh(lambda)              *  sin(lambda)             )/den;
k22 = (2*EI/L)  *lambda  *(sinh(lambda)*cosh(lambda) -  sin(lambda)* cos(lambda))/den;
k24 = (2*EI/L)  *lambda  *( sin(lambda)*cosh(lambda) - sinh(lambda)* cos(lambda))/den;

K2 = [  k11  k12 -k13  k14
        k12  k22 -k14  k24
       -k13 -k14  k11 -k12
        k14  k24 -k12  k22 ];

%% se verifica numéricamente que ambas matrices son iguales
valores.EI = 101265360;
valores.k  = 55241;         % escriba aquí un conjunto cualquiera de números
valores.L  = 5;

numK1 = real(double(subs(K1, valores)))
numK2 =      double(subs(K2, valores))

numK1 - numK2

%% Bye, bye!
