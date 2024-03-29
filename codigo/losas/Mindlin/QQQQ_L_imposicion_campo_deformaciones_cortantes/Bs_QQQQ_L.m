function Bbar_s = Bs_QQQQ_L(xi, eta, xe, ye, Nforma, dN_dxi, dN_deta, J_xi_eta)
%% Calcula la matriz de deformación sustitutiva por cortante Bs para el EF
%% de losa de Mindlin QQQQ-L
%
% Bbar_s = Bs_QQQQ_L(xi, eta, xe, ye, Nforma, dN_dxi, dN_deta, Je)
% 
% (xi, eta)        punto de integración de GL
% (xe, ye)         coordenadas de los nodos del EF

XI = 1; ETA = 2;

%% se define el número de nodos del EF
nno = 9;

%% Se definen los puntos de colocación
a = 1/sqrt(3);
nod = [ ...
% xi eta
   a   1    %  1
  -a   1    %  2
   a   0    %  3
  -a   0    %  4
   a  -1    %  5 
  -a  -1    %  6
   1   a    %  7
   1  -a    %  8          
   0   a    %  9
   0  -a    % 10          
  -1   a    % 11
  -1  -a ]; % 12
cxi  = nod(:,XI);
ceta = nod(:,ETA);
ngamma = length(cxi);     % número de puntos de colocación
 
Bbar_s = cell(ngamma,1);
J      = cell(ngamma,1);
for i = 1:ngamma
   %% Se evaluan las funciones de forma en los puntos de colocación
   N        = Nforma(cxi(i),  ceta(i));
   
   %% Se evalúan las derivadas de las funciones de forma en los puntos de colocación
   ddN_dxi  = dN_dxi (cxi(i), ceta(i));
   ddN_deta = dN_deta(cxi(i), ceta(i));

   %% Se utilizan las funciones de forma de w para el cálculo de la 
   %% transformación isoparamétrica   
   dx_dxi  = sum(ddN_dxi .*xe);   dy_dxi  = sum(ddN_dxi .*ye);
   dx_deta = sum(ddN_deta.*xe);   dy_deta = sum(ddN_deta.*ye);
   
   %% Se ensambla la matriz Jacobiana del elemento
   J{i} = [ dx_dxi   dy_dxi
            dx_deta  dy_deta ];

   %% Se calcula el determinante del Jacobiano
   det_Je = det(J{i});
   if det_Je <= 0
      error('El det_Je es negativo');
   end

   %% Se ensambla la matriz de deformación del elemento Bbar_s para el nodo i
   Bbar_s{i} = zeros(2,3*nno);
   for j = 1:nno
      % Se ensambla la matriz de deformación por cortante Bs
      % debo recalcular dN_dx y dN_dy para los puntos de colocación
      dN_dx_j = (+dy_deta*ddN_dxi(j) - dy_dxi*ddN_deta(j))/det_Je;
      dN_dy_j = (-dx_deta*ddN_dxi(j) + dx_dxi*ddN_deta(j))/det_Je;
      
      Bbar_s{i}(:, (3*j-2):(3*j)) = [ dN_dx_j  -N(j)  0     
                                      dN_dy_j   0    -N(j) ];                        
   end
end
Bhat_s = cell2mat(Bbar_s);               % eq 6.79
C = blkdiag(J{:});
       
% La matriz A_invP_T se dedujo con el programa APm1T_QQQQ_L_metodo1.m
A_invP_T = [ ...
  (eta*(3^(1/2)*xi + 1)*(eta + 1))/4,                                  0
                                   0,                                  0
 -(eta*(3^(1/2)*xi - 1)*(eta + 1))/4,                                  0
                                   0,                                  0
   -((eta^2 - 1)*(3^(1/2)*xi + 1))/2,                                  0
                                   0,                                  0
    ((eta^2 - 1)*(3^(1/2)*xi - 1))/2,                                  0
                                   0,                                  0
  (eta*(3^(1/2)*xi + 1)*(eta - 1))/4,                                  0
                                   0,                                  0
 -(eta*(3^(1/2)*xi - 1)*(eta - 1))/4,                                  0
                                   0,                                  0
                                   0,                                  0
                                   0,  (xi*(3^(1/2)*eta + 1)*(xi + 1))/4
                                   0,                                  0
                                   0, -(xi*(3^(1/2)*eta - 1)*(xi + 1))/4
                                   0,                                  0
                                   0,  -((xi^2 - 1)*(3^(1/2)*eta + 1))/2
                                   0,                                  0
                                   0,   ((xi^2 - 1)*(3^(1/2)*eta - 1))/2
                                   0,                                  0
                                   0,  (xi*(3^(1/2)*eta + 1)*(xi - 1))/4
                                   0,                                  0
                                   0, -(xi*(3^(1/2)*eta - 1)*(xi - 1))/4 ].';

Bbar_s = inv(J_xi_eta) * A_invP_T * C * Bhat_s;  % eq 6.80

return
