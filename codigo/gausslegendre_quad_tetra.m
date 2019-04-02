function [x_gl, w_gl] = gausslegendre_quad_tetra(n_gl)
%% Coordenadas y pesos de los puntos de integracion en las cuadraturas de 
%% Gauss-Legendre para el elementro tetraedrico

% Ver Onate pagina 312 (y Hughes pagina 174)
% Figura 7.21

switch n_gl
  case 1 % precision lineal
    x_gl = [1/4 1/4 1/4 1/4];
    w_gl = 1/6;    
    
  case 4 % precision cuadratica
    a = 0.58541020;
    b = 0.13819660;
    x_gl = [a b b b
            b a b b
            b b a b
            b b b a ];          
    w_gl = [1/24 1/24 1/24 1/24];
    
  case 5 % precision cubica
    g = -2/15;
    d = 3/40;
    x_gl = [1/4 1/4 1/4 1/4
            1/2 1/6 1/6 1/6   % OJO AQUI SE CORRIGIO UN ERROR DE ONATE
            1/6 1/2 1/6 1/6
            1/6 1/6 1/2 1/6
            1/6 1/6 1/6 1/2 ];
    w_gl = [g d d d d];
    
  otherwise 
    error('Orden de la cuadratura no soportada');
end

return
