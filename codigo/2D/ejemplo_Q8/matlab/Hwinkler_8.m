function He = Hwinkler_8(xnod, lado, kbalastro, espesor)
% Calcula la matriz de rigidez asociada a los apoyos del EF sobre una 
% cimentacion elastica. Sirve para un EF rectangular de 8 (serendipito) o 
% 9 (lagrangiano) nodos 
%
% SERENDIPITO 8          LAGRANGIANO 9
% xnod = [ x1e y1e       xnod = [ x1e y1e
%          x2e y2e                x2e y2e
%          ... ...                ... ...
%          x8e y8e ];             x9e y9e ];
%
% lado = 123, 345, 567, 781
%
% kbalastro = [ k1x k1y k2x k2y k3x k3y ]; % si apoya sobre el lado 123
%             [ k3x k3y k4x k4y k5x k5y ]; % si apoya sobre el lado 345
%             [ k5x k5y k6x k6y k7x k7y ]; % si apoya sobre el lado 567
%             [ k7x k7y k8x k8y k1x k1y ]; % si apoya sobre el lado 781

%% Se definen algunas constantes
X = 1; Y = 2;

%% Parametros de la cuadratura de Gauss-Legendre
n_gl = 5;                        % orden de la cuadratura de Gauss-Legendre
[x_gl, w_gl] = gausslegendre_quad(n_gl);

%% Se definen las funciones de forma unidimensionales y sus derivadas
NN = @(xi) [ ...
   xi*(xi-1)/2      % N1
   (1+xi)*(1-xi)    % N2
   xi*(1+xi)/2 ];   % N3

dNN_dxi = @(xi) [ ...
   xi - 1/2         % dN1_dxi
   -2*xi            % dN2_dxi
   xi + 1/2 ];      % dN3_dxi

%% Se definen los indices de los lados
switch lado
   case 123,  idx = [ 1 2 3 ];
   case 345,  idx = [ 3 4 5 ];
   case 567,  idx = [ 5 6 7 ];
   case 781,  idx = [ 7 8 1 ];      
   otherwise, error('Unicamente se permiten los lados 123, 345, 567 o 781');
end

%% Se calcula la integral
suma   = zeros(16);
N      = zeros(1,8);
dN_dxi = zeros(1,8);
for p = 1:n_gl
   N(idx)      = NN(x_gl(p));
   
   matN = [ N(1) 0    N(2) 0    N(3) 0    N(4) 0    N(5) 0    N(6) 0    N(7) 0    N(8) 0
            0    N(1) 0    N(2) 0    N(3) 0    N(4) 0    N(5) 0    N(6) 0    N(7) 0    N(8) ];
 
   dN_dxi(idx) = dNN_dxi(x_gl(p));
   
   dx_dxi = dN_dxi*xnod(:,X);
   dy_dxi = dN_dxi*xnod(:,Y);
   ds_dxi = hypot(dx_dxi, dy_dxi);

   suma = suma + matN'*matN*ds_dxi*w_gl(p);
end

% se calcula la matriz de rigidez asociada a los resortes sobre los lados
kb = zeros(16,1);
kb(reshape([2*idx-1; 2*idx],6,1)) = kbalastro(:);

He = espesor*suma*diag(kb);

%%
return;

