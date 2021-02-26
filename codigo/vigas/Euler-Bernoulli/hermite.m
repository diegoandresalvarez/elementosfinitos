%% Interpolacion polinomica de Hermite
clear, clc

%% Definicion de variables
syms xi u1 u2 u1p u2p
x1 = -1; x2 =  1;

%% Funciones de forma Lagrangianas
L1 = poly2sym(polyfit([-1 1], [1 0], 1), xi);
L2 = poly2sym(polyfit([-1 1], [0 1], 1), xi);

%% Se calculan los polinomios Pi y Qi
P1 = (1 - 2*subs(diff(L1, xi), xi, x1)*(xi - x1))*L1^2;
P2 = (1 - 2*subs(diff(L2, xi), xi, x2)*(xi - x2))*L2^2;
Q1 = (xi - x1)*L1^2; 
Q2 = (xi - x2)*L2^2;

%% Se calcula el polinomio interpolador
u = simplify(P1*u1 + P2*u2 + Q1*u1p + Q2*u2p);

%% Se calculan las funciones de forma
N1  = expand(subs(u, {u1, u1p, u2, u2p}, {1, 0, 0, 0})); % = P1
N1b = expand(subs(u, {u1, u1p, u2, u2p}, {0, 1, 0, 0})); % = Q1
N2  = expand(subs(u, {u1, u1p, u2, u2p}, {0, 0, 1, 0})); % = P2
N2b = expand(subs(u, {u1, u1p, u2, u2p}, {0, 0, 0, 1})); % = Q2

%% Se muestran los resultados
disp('N_1(xi)  = '); pretty(N1);    
disp('Nb_1(xi) = '); pretty(N1b);
disp('N_2(xi)  = '); pretty(N2);    
disp('Nb_2(xi) = '); pretty(N2b);

%% Se dibujan las funciones de forma Hermitianas
figure
hold on
fplot(N1,  [-1,1], 'LineWidth', 2);
fplot(N1b, [-1,1], 'LineWidth', 2);
fplot(N2,  [-1,1], 'LineWidth', 2);
fplot(N2b, [-1,1], 'LineWidth', 2);
xlabel('\xi', 'FontSize', 20)
legend('$N_1(\xi)$', '$\bar{N}_1(\xi)$', '$N_2(\xi)$', '$\bar{N}_2(\xi)$', ...
       'Location', 'Best',     ...
       'Interpreter', 'Latex', ...
       'FontSize', 20);
title('Funciones de forma de la viga de Euler-Bernoulli de dos nodos')
axis equal
grid on

%% bye bye
