function [Bf, Bbar_c, det_J_xi_eta] = Bf_Bc_QQQQ_L(xi, eta, xe, ye)
%% Calcula la matriz de deformacion por flexion Bf y la matriz de deformacion
%% sustitutiva por cortante Bc para el elemento de losa de RM QQQQ-L
%
% [Bf, Bbar_c, det_J_xi_eta] = Bf_Bc_QQQQ_L(xi, eta, xe, ye)
% 
% (xi, eta)        punto de integracion de GL
% (xe, ye)         coordenadas de los nodos del EF

X = 1; Y = 2;

%% Se definen las funciones de forma
c9_funciones_forma_lagrangiano_9_nodos  % Nforma, dN_dxi, dN_deta
nno = 9;                                % numero de nodos del elemento finito

%% Se calcula el Jacobiano y su inversa para el punto de integracion (xi,eta)
ddN_dxi = dN_dxi(xi,eta);      ddN_deta = dN_deta(xi,eta);
dx_dxi   = sum(ddN_dxi .*xe);   dy_dxi    = sum(ddN_dxi .*ye);
dx_deta  = sum(ddN_deta.*xe);   dy_deta   = sum(ddN_deta.*ye);

% Se ensambla la matriz Jacobiana del elemento
J_xi_eta = [ dx_dxi   dy_dxi
             dx_deta  dy_deta ];
       
inv_J_xi_eta = inv(J_xi_eta);

% Se calcula el determinante del Jacobiano
det_J_xi_eta =  det(J_xi_eta);
if det_J_xi_eta <= 0
   error('El det_J_xi_eta es negativo');
end

%% Se calcula la matriz Bf
Bf = zeros(3,3*nno);
dN_dx = zeros(1,nno);
dN_dy = zeros(1,nno);
for i = 1:nno
    % Se ensambla la matriz de deformacion por flexion Bf
    dN_dx(i) = inv_J_xi_eta(1,1)*ddN_dxi(i) + inv_J_xi_eta(1,2)*ddN_deta(i);
    dN_dy(i) = inv_J_xi_eta(2,1)*ddN_dxi(i) + inv_J_xi_eta(2,2)*ddN_deta(i);
    Bf(:,[3*i-2 3*i-1 3*i]) = [ 0  -dN_dx(i)          0    % aqui se ensambla
                                0          0  -dN_dy(i)    % y asigna la matriz
                                0  -dN_dy(i)  -dN_dx(i) ]; % Bf_i
end;

%% Se definen los puntos de colocacion
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
cx = nod(:,X);
cy = nod(:,Y);
npc = length(cx);     % numero de puntos de colocacion
 
Bbar_c = cell(npc,1);
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

   %% Se ensambla la matriz de deformacion del elemento Bbar_c para el nodo i
   Bbar_c{i} = zeros(2,3*nno);
   for j = 1:nno
      % Se ensambla la matriz de deformacion por cortante Bc
      % debo recalcular dN_dx y dN_dy para los puntos de colocacion
      dN_dx_j = (+dy_deta*ddN_dxi(j) - dy_dxi*ddN_deta(j))/det_Je;
      dN_dy_j = (-dx_deta*ddN_dxi(j) + dx_dxi*ddN_deta(j))/det_Je;
      
      Bbar_c{i}(:,[(3*j-2):(3*j)]) = [ dN_dx_j  -N(j)  0     
                                       dN_dy_j   0    -N(j) ];                        
   end;
end
Bhat_c = cell2mat(Bbar_c);               % eq 6.79
C = blkdiag(J{:});
       
% La matriz A_invP_T se dedujo con el programa APm1T_QQQQ_L_metodo1.m
A_invP_T = [ ...
        (3^(1/2)/4)*eta*(xi + 3^(1/2)/3)*(eta + 1),                                                0
                                                 0,                                                0
       (-3^(1/2)/4)*eta*(xi - 3^(1/2)/3)*(eta + 1),                                                0
                                                 0,                                                0
 (-3^(1/2)/2)*(eta - 1)*(eta + 1)*(xi + 3^(1/2)/3),                                                0
                                                 0,                                                0
  (3^(1/2)/2)*(eta - 1)*(eta + 1)*(xi - 3^(1/2)/3),                                                0
                                                 0,                                                0
        (3^(1/2)/4)*eta*(xi + 3^(1/2)/3)*(eta - 1),                                                0
                                                 0,                                                0
       (-3^(1/2)/4)*eta*(xi - 3^(1/2)/3)*(eta - 1),                                                0
                                                 0,                                                0
                                                 0,                                                0
                                                 0,        (3^(1/2)/4)*xi*(eta + 3^(1/2)/3)*(xi + 1)
                                                 0,                                                0
                                                 0,       (-3^(1/2)/4)*xi*(eta - 3^(1/2)/3)*(xi + 1)
                                                 0,                                                0
                                                 0, (-3^(1/2)/2)*(xi - 1)*(xi + 1)*(eta + 3^(1/2)/3)
                                                 0,                                                0
                                                 0,  (3^(1/2)/2)*(xi - 1)*(xi + 1)*(eta - 3^(1/2)/3)
                                                 0,                                                0
                                                 0,        (3^(1/2)/4)*xi*(eta + 3^(1/2)/3)*(xi - 1)
                                                 0,                                                0
                                                 0,       (-3^(1/2)/4)*xi*(eta - 3^(1/2)/3)*(xi - 1) ]';

Bbar_c = inv_J_xi_eta * A_invP_T * C * Bhat_c;  % eq 6.80
    
return
