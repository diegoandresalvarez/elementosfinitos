function [B, Je, det_J] = B_lam_deg(e, xi, eta, zeta)
%% Calcula la matriz B para el elemento de lamina degenerada en el punto 
%% de coordenadas (xi,eta,zeta)
%
% [B, Je, det_J] = B_lam_deg(e, xi, eta, zeta)
% 
% e                    numero del EF
% (xi, eta, zeta)      punto de integracion de GL

global Nforma dN_dxi dN_deta vg1 vg2 vg3 t LaG xnod

LaGe = LaG(e,:);     % nodos del EF
vg1e = vg1(LaGe,:);  % vectores que definen el sistema de coord nodales
vg2e = vg2(LaGe,:);  % vgi = [nnoef x 3], para i = 1,2,3
vg3e = vg3(LaGe,:);
r0e  = xnod(LaGe,:); % [xe, ye, ze] = [nnoef x 3] coord de los nodos del EF
te   = t(LaGe);      % [nnoef x 1] espesor en cada uno de los nodos del EF

%% Se calculan las funciones de forma y sus derivadas
NN       = Nforma(xi,eta);
dNN_dxi  = dN_dxi(xi,eta);
dNN_deta = dN_deta(xi,eta);

%% Se calcula el numero de nodos del EF
nnoef = length(LaGe);

%% Se calcula zb
zb = zeta*te/2; %(zeta - zeta0)*t/2;

%% Se calcula dr_dxi, dr_deta y dr_dzeta, las cuales son las filas del 
%% Jacobiano
dr_dxi   = zeros(1,3);
dr_deta  = zeros(1,3);
dr_dzeta = zeros(1,3);

for i = 1:nnoef
   ri = r0e(i,:) + zb(i)*vg3e(i,:);
   dr_dxi   =  dr_dxi   + dNN_dxi (i)*ri;
   dr_deta  =  dr_deta  + dNN_deta(i)*ri;
   dr_dzeta =  dr_dzeta + NN(i)*vg3e(i,:)*te(i)/2;
end

%% Se calcula el Jacobiano para el punto de integracion (xi,eta)
% Se ensambla la matriz Jacobiana del elemento
Je = [ dr_dxi; dr_deta; dr_dzeta ];

% Se calcula el determinante del Jacobiano y se verifica que sea positivo
det_J = det(Je);
if det_J >= 0
   error('El det_J es positivo o cero');
end

%% Se calcula la inv(Je)
invJe = inv(Je);

invJ_11 = invJe(1,1);   invJ_12 = invJe(1,2);   invJ_13 = invJe(1,3);
invJ_21 = invJe(2,1);   invJ_22 = invJe(2,2);   invJ_23 = invJe(2,3);
invJ_31 = invJe(3,1);   invJ_32 = invJe(3,2);   invJ_33 = invJe(3,3);

%% Se calcula la matriz B (ver: c12_matriz_B.m)
B = zeros(6,5*nnoef);
for i = 1:nnoef
   zi       = zb(i);
   ti       = te(i);
   Ni       = NN(i);
   dNi_dxi  = dNN_dxi(i);
   dNi_deta = dNN_deta(i);

   % Se calculan los terminos de Ci
   Ci = [vg1e(i,:)' vg2e(i,:)'];

   Ni_1 = dNi_dxi*invJ_11 + dNi_deta*invJ_12;
   Ni_2 = dNi_dxi*invJ_21 + dNi_deta*invJ_22;
   Ni_3 = dNi_dxi*invJ_31 + dNi_deta*invJ_32;
   Gi_1 = -(ti*invJ_13*Ni/2 + zi*Ni_1)*Ci;
   Gi_2 = -(ti*invJ_23*Ni/2 + zi*Ni_2)*Ci;
   Gi_3 = -(ti*invJ_33*Ni/2 + zi*Ni_3)*Ci;

   B(:,(5*i-4):(5*i)) = [ ...
     Ni_1,    0,    0,               Gi_1(1,1),               Gi_1(1,2)
        0, Ni_2,    0,               Gi_2(2,1),               Gi_2(2,2)
        0,    0, Ni_3,               Gi_3(3,1),               Gi_3(3,2)
     Ni_2, Ni_1,    0,   Gi_2(1,1) + Gi_1(2,1),   Gi_2(1,2) + Gi_1(2,2)
     Ni_3,    0, Ni_1,   Gi_3(1,1) + Gi_1(3,1),   Gi_3(1,2) + Gi_1(3,2)     
        0, Ni_3, Ni_2,   Gi_3(2,1) + Gi_2(3,1),   Gi_3(2,2) + Gi_2(3,2) ];   
end 

return
