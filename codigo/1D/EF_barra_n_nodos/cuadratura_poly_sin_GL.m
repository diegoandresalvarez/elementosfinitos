clear, clc, close all

% Integracion de:
%
% f(x) = 0.2 + 25*x - 200*x.^2 +675*x.^3 - 900*x.^4 + 400*x.^5
%
% entre 0 y 0.8 usando cuadraturas de Gauss-Legendre

a   = 0;                  % Limites de integracion
b   = 0.8;        
f   = @(x) 0.2 + 25*x - 200*x.^2 +675*x.^3 - 900*x.^4 + 400*x.^5; % la funcion
sol = 3076/1875;
err = zeros(10,1);        % Separo la memoria
for m = 1:10;             % Vario el numero de puntos de la cuadratura
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
% f(x) = sin(x)
%
% entre 0 y pi/2 usando cuadraturas de Gauss-Legendre

a   = 0;
b   = pi/2;
f   = @(x) sin(x);
sol = 1;
err = zeros(10,1);
for m = 1:10;
   [x,c] = gausslegendre_quad(m);
   err(m) = abs(((b-a)/2)*sum(c.*f((b+a)/2 + (b-a)*x/2)) - sol);
end
figure
semilogy(err); % se grafica en escala logaritmica para apreciar el error
xlabel('Numero de puntos en la cuadratura');
ylabel('Error');
title('Cuadratura de Gauss Legendre');
grid
