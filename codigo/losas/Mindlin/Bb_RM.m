function [Bb, det_Je, Je] = Bb_RM(xi, eta, xe, ye, dN_dxi, dN_deta)
%% Calcula la matriz de deformaciones por flexion Bb
%
% [Bb, det_Je, Je] = Bb_RM(xi, eta, xe, ye, dN_dxi, dN_deta)
%
% (xi, eta)        punto de integracion de GL
% (xe, ye)         coordenadas de los nodos del EF
% dN_dxi, dN_deta  function handle derivadas de las funciones de forma del EF

%% Se evaluan las derivadas de las funciones de forma en los puntos
%% de integracion de Gauss-Legendre
ddN_dxi  = dN_dxi (xi, eta);
ddN_deta = dN_deta(xi, eta);

%% Se utilizan las funciones de forma de w para el calculo de la 
%% transformacion isoparametrica
dx_dxi  = sum(ddN_dxi  .* xe);   dy_dxi  = sum(ddN_dxi  .* ye);
dx_deta = sum(ddN_deta .* xe);   dy_deta = sum(ddN_deta .* ye);

%% Se ensambla la matriz Jacobiana del elemento
Je = [ dx_dxi   dy_dxi
       dx_deta  dy_deta ];

%% Se calcula el determinante del Jacobiano
det_Je = det(Je);
if det_Je <= 0
   error('El det_Je es negativo');
end

%% Se ensambla la matriz de deformacion del elemento Bb
nno = length(xe);
Bb = zeros(3,3*nno);
for i = 1:nno   
   dNi_dx = (+dy_deta*ddN_dxi(i) - dy_dxi*ddN_deta(i))/det_Je;
   dNi_dy = (-dx_deta*ddN_dxi(i) + dx_dxi*ddN_deta(i))/det_Je;

   Bb(:,(3*i-2):(3*i)) = [ 0 -dNi_dx       0    % se ensambla y
                           0       0 -dNi_dy    % asigna la matriz
                           0 -dNi_dy -dNi_dx ]; % Bb_i
end

return
