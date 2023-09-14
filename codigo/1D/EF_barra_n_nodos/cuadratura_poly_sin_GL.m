clear, clc, close all

% Integración de:
f   = @(x) 0.2 + 25*x - 200*x.^2 + 675*x.^3 - 900*x.^4 + 400*x.^5;
% entre 0 y 0.8 usando cuadraturas de Gauss-Legendre

a   = 0;                  % límites de integración
b   = 0.8;        
sol = 3076/1875;          % solución exacta
err = zeros(10,1);        % se reserva la memoria
for m = 1:10              % se varía el número de puntos de la cuadratura
   [xi,w] = gausslegendre_quad(m);  % cálculo de las raíces y los pesos
   err(m) = abs(((b-a)/2)*sum(w.*f((b+a)/2 + (b-a)*xi/2)) - sol);
end
figure                    % creo un lienzo
plot(err)                 % grafico el error
xlabel('Número de puntos en la cuadratura');
ylabel('Error absoluto');
title('Cuadratura de Gauss-Legendre');
grid minor                % pongo la rejilla


% Integración de:
f   = @(x) sin(x);
% entre 0 y pi/2 usando cuadraturas de Gauss-Legendre

a   = 0;                  % límites de integración
b   = pi/2;
sol = 1;                  % solución exacta
err = zeros(10,1);        % se reserva la memoria
for m = 1:10              % se varía el número de puntos de la cuadratura
   [xi,w] = gausslegendre_quad(m);  % cálculo de las raíces y los pesos
   err(m) = abs(((b-a)/2)*sum(w.*f((b+a)/2 + (b-a)*xi/2)) - sol);
end
figure
semilogy(err); % escala logarítmica para apreciar mejor el error
xlabel('Número de puntos en la cuadratura');
ylabel('Error absoluto');
title('Cuadratura de Gauss-Legendre');
grid minor

