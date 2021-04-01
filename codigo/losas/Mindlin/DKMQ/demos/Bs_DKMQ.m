function Bbar_s = Bs_QQQQ_L(xi, eta, xe, ye, J_xi_eta)
%% Calcula la matriz de deformacion sustitutiva por cortante Bs para el EF
%% de losa de Mindlin QQQQ-L
%
% Bbar_s = Bs_QQQQ_L(xi, eta, xe, ye, Nforma, dN_dxi, dN_deta, Je)
% 
% (xi, eta)        punto de integracion de GL
% (xe, ye)         coordenadas de los nodos del EF

X = 1; Y = 2;

% Funciones de forma serendipitas del elemento rectangular de 4 nodos:
Nforma = @(xi,eta) [ ...
     ((eta - 1)*(xi - 1))/4
    -((eta - 1)*(xi + 1))/4
     ((eta + 1)*(xi + 1))/4
    -((eta + 1)*(xi - 1))/4 ];

% Derivadas con respecto a xi:
dN_dxi = @(xi,eta) [ ...
       eta/4 - 1/4
       1/4 - eta/4
       eta/4 + 1/4
     - eta/4 - 1/4 ];

% Derivadas con respecto a eta:
dN_deta = @(xi,eta) [ ...
       xi/4 - 1/4
     - xi/4 - 1/4
       xi/4 + 1/4
       1/4 - xi/4 ];

%% se define el numero de nodos del EF
nno = 4;

%% Se definen los puntos de colocacion
nod = [ ...
% xi eta
   0  -1   %  5
   1   0   %  6
   0   1   %  7
  -1   0 ]; % 8
cx = nod(:,X);
cy = nod(:,Y);
npc = length(cx);     % numero de puntos de colocacion
 
Bbar_s = cell(npc,1);
J      = cell(npc,1);
for i = 1:npc
   %% Se evaluan las funciones de forma en los puntos de colocacion   
   N        = Nforma(cx(i),  cy(i));
   
   %% Se evaluan las derivadas de las funciones de forma en los puntos de colocacion   
   ddN_dxi  = dN_dxi (cx(i), cy(i));
   ddN_deta = dN_deta(cx(i), cy(i));

   %% Se utilizan las funciones de forma de w para el calculo de la 
   %% transformacion isoparametrica   
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

   %% Se ensambla la matriz de deformacion del elemento Bbar_s para el nodo i
   Bbar_s{i} = zeros(2,3*nno);
   for j = 1:nno
      % Se ensambla la matriz de deformacion por cortante Bs
      % debo recalcular dN_dx y dN_dy para los puntos de colocacion
      dN_dx_j = (+dy_deta*ddN_dxi(j) - dy_dxi*ddN_deta(j))/det_Je;
      dN_dy_j = (-dx_deta*ddN_dxi(j) + dx_dxi*ddN_deta(j))/det_Je;
      
      Bbar_s{i}(:, (3*j-2):(3*j)) = [ dN_dx_j  -N(j)  0     
                                      dN_dy_j   0    -N(j) ];                        
   end
end
Bhat_s = cell2mat(Bbar_s);               % eq 6.79
C = blkdiag(J{:});
       
% La matriz A_invP_T se dedujo con el programa APm1T_QQQQ_L_metodo1.m
A_invP_T =[ 1/2 - eta/2, 0, 0,          0, eta/2 + 1/2, 0, 0,          0
                      0, 0, 0, xi/2 + 1/2,           0, 0, 0, 1/2 - xi/2 ].';

J_xi_eta = [ j11 j12; j21 j22 ];

Bbar_s = inv(J_xi_eta) * A_invP_T * C * Bhat_s;  % eq 6.80
    
return
