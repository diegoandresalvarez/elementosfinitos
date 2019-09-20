clear, clc, close all

% Integracion de:
%
f   = @(x) 0.2 + 25*x - 200*x.^2 +675*x.^3 - 900*x.^4 + 400*x.^5;
%
% entre 0 y 0.8 usando cuadraturas de Gauss-Legendre

a   = 0;                  % limites de integracion
b   = 0.8;        
sol = 3076/1875;          % solucion exacta
err = zeros(10,1);        % separo la memoria
for m = 1:10              % vario el numero de puntos de la cuadratura
   [x,c] = gausslegendre_quad(m);  % calculo w y c de la cuadratura
   err(m) = abs(((b-a)/2)*sum(c.*f((b+a)/2 + (b-a)*x/2)) - sol);
end
figure                    % creo un lienzo
plot(err)                 % grafico el error
xlabel('Numero de puntos en la cuadratura');
ylabel('Error');
title('Cuadratura de Gauss Legendre');
grid                      % pongo la rejilla


% Integracion de:
%
f   = @(x) sin(x);
%
% entre 0 y pi/2 usando cuadraturas de Gauss-Legendre

a   = 0;                  % limites de integracion
b   = pi/2;
sol = 1;                  % solucion exacta
err = zeros(10,1);        % separo la memoria
for m = 1:10
   [x,c] = gausslegendre_quad(m);
   err(m) = abs(((b-a)/2)*sum(c.*f((b+a)/2 + (b-a)*x/2)) - sol);
end
figure
semilogy(err); % escala logaritmica para apreciar mejor el error
xlabel('Numero de puntos en la cuadratura');
ylabel('Error');
title('Cuadratura de Gauss Legendre');
grid
