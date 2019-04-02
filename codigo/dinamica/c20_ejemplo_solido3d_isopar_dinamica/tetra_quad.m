function [L2L3L4,W] = tetra_quad(n)
% Cuadratura de Gauss-Legendre para un tetraedro

switch n
   case 1 % Precision lineal
      L2L3L4 = [0.25 0.25 0.25];
      W = 1;
   case 4 % Precision cuadratica
      a = 0.58541020; b = 0.13819660;
      L2L3L4 = [ a b b
                 b a b
                 b b a
                 b b b ];
      W = [0.25; 0.25; 0.25; 0.25];
   case 5 % Precision cubica
      L2L3L4 = [0.25 0.25 0.25
                1/3  1/6  1/6
                1/6  1/3  1/6
                1/6  1/6  1/3
                1/6  1/6  1/6 ];      
      W = [-0.8; 0.45; 0.45; 0.45; 0.45];
   otherwise
      error('Numero de puntos de la cuadratura no soportada');         
end
W = W/6;

return; % bye, bye
