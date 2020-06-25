function [Bf, Bm, Bbar_c, det_J_xi_eta] = Bf_Bm_Bc_QLLL(xi, eta, xe, ye, Te)
%% Calcula las matrices Bf, Bc y Bm para el elemento de lamina de RM QLLL
%
% [Bf, Bm, Bbar_c, det_J_xi_eta] = Bf_Bm_Bc_QLLL(xi, eta, xe, ye, Te)
% 
% (xi, eta)        punto de integracion de GL
% (xe, ye)         coordenadas de los nodos del EF

X = 1; Y = 2;

%% Se definen las funciones de forma
c9_funciones_forma_rect_4_nodos  % Nforma, dN_dxi, dN_deta
nno = 4;                         % numero de nodos del elemento finito

%% Se calcula el Jacobiano y su inversa para el punto de integracion (xi,eta)
ddN_dxi = dN_dxi(xi,eta);      ddN_deta = dN_deta(xi,eta);
dx_dxi  = sum(ddN_dxi .*xe);   dy_dxi   = sum(ddN_dxi .*ye);
dx_deta = sum(ddN_deta.*xe);   dy_deta  = sum(ddN_deta.*ye);

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
Bf = zeros(3,5*nno);
Bm = zeros(3,5*nno);
dN_dx = zeros(1,nno);
dN_dy = zeros(1,nno);
for i = 1:nno
    % Se ensambla la matriz de deformacion por flexion Bf
    dN_dx(i) = inv_J_xi_eta(1,1)*ddN_dxi(i) + inv_J_xi_eta(1,2)*ddN_deta(i);
    dN_dy(i) = inv_J_xi_eta(2,1)*ddN_dxi(i) + inv_J_xi_eta(2,2)*ddN_deta(i);
    Bf(:,(5*i-4):(5*i)) = [ 0 0 0  -dN_dx(i)          0    % aqui se ensambla
                            0 0 0          0  -dN_dy(i)    % y asigna la matriz
                            0 0 0  -dN_dy(i)  -dN_dx(i) ]; % Bf_i
                             
    Bm(:,(5*i-4):(5*i)) = [ dN_dx(i)         0 0 0 0       % aqui se ensambla
                                   0  dN_dy(i) 0 0 0       % y asigna la matriz
                            dN_dy(i)  dN_dx(i) 0 0 0];     % Bm_i
end;

%% Se definen los puntos de colocacion
cx = [  0;  1;  0; -1 ];
cy = [ -1;  0;  1;  0 ];

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
      
      Bbar_c{i}(:,(3*j-2):(3*j)) = [ dN_dx_j  -N(j)  0     
                                     dN_dy_j   0    -N(j) ];                        
   end;
end
Bhat_c = cell2mat(Bbar_c);               % eq 6.79
C = blkdiag(J{:});
       
% La matriz A_invP_T se dedujo con el programa APm1T_QLLL_metodo1.m
A_invP_T = [ ...
 1/2-eta/2, 0, 0,          0, eta/2+1/2, 0, 0,        0
         0, 0, 0, xi/2+1/2,           0, 0, 0, 1/2-xi/2 ];
      
Bbar_c = inv_J_xi_eta * A_invP_T * C * Bhat_c;  % eq 6.80
    
% Se expande la matriz con ceros
tmp = zeros(2,5*nno);
for j = 1:nno
   tmp(:,(5*j-4):(5*j)) = [ zeros(2,2) Bbar_c(:,(3*j-2):(3*j)) ]; 
end;
Bbar_c = tmp;

%% Se aplican las matrices de transformacion
Bf     = Bf*Te;
Bm     = Bm*Te;
Bbar_c = Bbar_c*Te;

return
