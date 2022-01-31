clear, clc, close all
%% -------------------------------------------------------------------------
%% Funciones de Forma Lagrangianas de dos nodos

% Calculo las funciones de forma
syms xi
N1 = calc_N([-1 1], [1 0], xi);
N2 = calc_N([-1 1], [0 1], xi);

% Imprimo las funciones de forma
fprintf('\n\nFunciones de Forma Lagrangianas de DOS nodos:\n')
fprintf('\n\nN1(xi) = \n'), pretty(N1)
fprintf('\n\nN2(xi) = \n'), pretty(N2)

% Grafico las funciones de forma
figure                 % Creo un lienzo
grid on                % creo la rejilla
hold on;               % Para que no se sobreescriban los graficos
h1 = fplot(N1, [-1 1], 'Color', 'r', 'LineWidth', 2);
h2 = fplot(N2, [-1 1], 'Color', 'b', 'LineWidth', 2);
title('Funciones de Forma Lagrangianas de DOS nodos')
xlabel('\xi'); 
ylabel('N_i(\xi)');
plot([-1 1],[0 0], 'ko', [-1 1],[1 1], 'ko') % grafico los nodos
legend([h1, h2], 'N1(\xi)','N2(\xi)','Location','Best');

%% -------------------------------------------------------------------------
%% Funciones de Forma Lagrangianas de tres nodos

% Calculo las funciones de forma
syms xi
N1 = calc_N([-1 0 1], [1 0 0], xi); % = xi*(xi-1)/2
N2 = calc_N([-1 0 1], [0 1 0], xi); % = (1+xi)*(1-xi)
N3 = calc_N([-1 0 1], [0 0 1], xi); % = xi*(xi+1)/2

% Imprimo las funciones de forma
fprintf('\n\nFunciones de Forma Lagrangianas de TRES nodos:\n')
fprintf('\n\nN1(xi) = \n'), pretty(N1)
fprintf('\n\nN2(xi) = \n'), pretty(N2)
fprintf('\n\nN3(xi) = \n'), pretty(N3)

% Grafico las funciones de forma
figure                 % Creo un lienzo
grid on                % creo la rejilla
hold on;               % Para que no se sobreescriban los graficos
h1 = fplot(N1, [-1 1], 'Color', 'r', 'LineWidth', 2);
h2 = fplot(N2, [-1 1], 'Color', 'b', 'LineWidth', 2);
h3 = fplot(N3, [-1 1], 'Color', 'c', 'LineWidth', 2);
title('Funciones de Forma Lagrangianas de TRES nodos')
xlabel('\xi'); 
ylabel('N_i(\xi)');
plot([-1 0 1],[0 0 0], 'ko', [-1 0 1],[1 1 1], 'ko') % grafico los nodos
legend([h1, h2, h3], 'N1(\xi)','N2(\xi)','N3(\xi)','Location','Best');

%% -------------------------------------------------------------------------
%% Funciones de Forma Lagrangianas de cuatro nodos

% Calculo las funciones de forma
syms xi
N1 = calc_N([-1 -1/3 1/3 1], [1 0 0 0], xi);
N2 = calc_N([-1 -1/3 1/3 1], [0 1 0 0], xi);
N3 = calc_N([-1 -1/3 1/3 1], [0 0 1 0], xi);
N4 = calc_N([-1 -1/3 1/3 1], [0 0 0 1], xi);

% Imprimo las funciones de forma
fprintf('\n\nFunciones de Forma Lagrangianas de CUATRO nodos:\n')
fprintf('\n\nN1(xi) = \n'), pretty(N1);  fprintf('\n\nN2(xi) = \n'), pretty(N2)
fprintf('\n\nN3(xi) = \n'), pretty(N3);  fprintf('\n\nN4(xi) = \n'), pretty(N4)

% Grafico las funciones de forma
figure                 % Creo un lienzo
grid on                % creo la rejilla
hold on;               % Para que no se sobreescriban los graficos
h1 = fplot(N1, [-1 1], 'Color', 'r', 'LineWidth', 2);
h2 = fplot(N2, [-1 1], 'Color', 'b', 'LineWidth', 2);
h3 = fplot(N3, [-1 1], 'Color', 'c', 'LineWidth', 2);
h4 = fplot(N4, [-1 1], 'Color', 'm', 'LineWidth', 2);
title('Funciones de Forma Lagrangianas de CUATRO nodos')
xlabel('\xi'); 
ylabel('N_i(\xi)');
plot([-1 -1/3 1/3 1],[0 0 0 0], 'ko', [-1 -1/3 1/3 1],[1 1 1 1], 'ko')
axis([-1 1 -0.4 1.2])
legend([h1 h2 h3 h4], 'N1(\xi)','N2(\xi)','N3(\xi)','N4(\xi)','Location','Best');

%% -------------------------------------------------------------------------
%% Funciones de Forma Lagrangianas de cinco nodos

% Calculo las funciones de forma
syms xi
N = cell(5,1);
N{1} = calc_N([-1 -1/2 0 1/2 1], [1 0 0 0 0], xi);
N{2} = calc_N([-1 -1/2 0 1/2 1], [0 1 0 0 0], xi);
N{3} = calc_N([-1 -1/2 0 1/2 1], [0 0 1 0 0], xi);
N{4} = calc_N([-1 -1/2 0 1/2 1], [0 0 0 1 0], xi);
N{5} = calc_N([-1 -1/2 0 1/2 1], [0 0 0 0 1], xi);

% Imprimo las funciones de forma
fprintf('\n\nFunciones de Forma Lagrangianas de CINCO nodos:\n')
for i = 1:5
    fprintf('\n\nN%d(xi) = \n', i); pretty(N{i});
end

% Grafico las funciones de forma
figure                 % Creo un lienzo
grid on                % creo la rejilla
hold on;               % Para que no se sobreescriban los graficos
h = zeros(5,1);
color = ['r', 'b', 'c', 'm', 'k'];
for i = 1:5
    h(i) = fplot(N{i}, [-1 1], 'Color', color(i), 'LineWidth', 2);
end
title('Funciones de Forma Lagrangianas de CINCO nodos')
xlabel('\xi'); 
ylabel('N_i(\xi)');
plot([-1 -1/2 0 1/2 1],[0 0 0 0 0], 'ko', [-1 -1/2 0 1/2 1],[1 1 1 1 1], 'ko')
axis([-1 1 -0.6 1.2])
legend(h, 'N1(\xi)','N2(\xi)','N3(\xi)','N4(\xi)','N5(\xi)', 'Location','Best');

%% -------------------------------------------------------------------------
%% Calcular correctamente los polinomios de las funciones de forma 1D
function N = calc_N(xp, yp, var)
    % se ve verifican los tamanios de los vectores xp y yp
    nx = length(xp);
    ny = length(yp);
    assert(nx == ny, 'Los vectores xp y yp deben tener el mismo tamanio');

    % se calculan los coeficientes de los polinomios
    c = polyfit(xp, yp, nx-1);
    
    % se eliminan los errores en la aproximacion numerica, haciendo los
    % coeficientes demasiado pequenios igual a cero
    c(abs(c) < 1e-10) = 0;
    
    % con los coeficientes corregidos se calculan las funciones de forma
    N = poly2sym(c, var);
end
