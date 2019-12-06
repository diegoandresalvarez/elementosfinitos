%% PORTICO 3D
% Topologia tomada de:
% Van-Long Hoang, Hung Nguyen Dang, Jean-Pierre Jaspart and Jean-François 
% Demonceau (2015). An overview of the plastic-hinge analysis of 3D steel 
% frames. Asia Pacific Journal on Computational Engineering. Vol. 2. No. 4.
% https://doi.org/10.1186/s40540-015-0016-9

% Considere el pórtico de seis pisos mostrado en la figura. El módulo de 
% elasticidad en todos los miembros es 206 x 10^3 N/mm^2, 
% nu(Poisson) = 0.30. La carga en cada piso es uniforme e igual a 
% 4.8 kN/m^2; las cargas de viento se simulan como cargas puntuales de 
% 26.7 kN en la dirección Y en cada unión viga-columna.

clear, clc, close all

%% constantes
X = 1; Y = 2; Z = 3; 
U = 1; V = 2; W = 3; THX = 4; THY = 5; THZ = 6;

%% Unidades en kN y mm

% Cargar:
% * xnod: coordenadas de cada nodo [x, y, z]
% * LaG: local a global: matriz que relaciona nodos locales y globales
% * V3_V1
% * tipo_barra
% * propiedades (A(mm2)	Iz(mm4)	Iy(mm4)	J(mm4))
load portico_3d;

A  = propiedades(tipo_barra, 1); % mm2
Iz = propiedades(tipo_barra, 2); % mm4
Iy = propiedades(tipo_barra, 3); % mm4
J  = propiedades(tipo_barra, 4); % mm4

E  = 206;          % modulo de elasticidad (kN/mm^2) 
nu = 0.30;         % modulo de Poisson
G  = E/(2*(1+nu)); % modulo de rigidez

nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
nbar = size(LaG,1);  % numero de EFs (numero de filas de LaG)
ngdl = 6*nno;        % numero de grados de libertad (tres por nodo)

%% gdl: grados de libertad
% fila = nodo
% col1 = gdl desplazamiento en direccion x
% col2 = gdl desplazamiento en direccion y
% col3 = gdl desplazamiento en direccion z
% col4 = gdl rotacion alrededor de eje x
% col5 = gdl rotacion alrededor de eje x
% col6 = gdl rotacion alrededor de eje x
gdl = reshape(1:ngdl,6,nno)';

%% cargas aplicadas
nodos_cargados_viento = find(xnod(:,Y) == 0)';
dofs_cargados_viento  = 6*nodos_cargados_viento - 4;  % cargas en la dir. Y

f = zeros(ngdl, 1);
f(dofs_cargados_viento) = 26.7; %kN

%% Se dibuja la estructura junto con su numeracion
figure; 
hold on;
box on
view(3)
for e = 1:nbar
   line(xnod(LaG(e,:),X), xnod(LaG(e,:),Y) , xnod(LaG(e,:),Z));
   
   % Calculo la posicion del centro de gravedad del triangulo
   cgx = (xnod(LaG(e,1),X) + xnod(LaG(e,2),X))/2;
   cgy = (xnod(LaG(e,1),Y) + xnod(LaG(e,2),Y))/2;   
   cgz = (xnod(LaG(e,1),Z) + xnod(LaG(e,2),Z))/2;      
   h = text(cgx, cgy, cgz, num2str(e)); set(h,'Color', [1 0 0]);
end
axis equal
grid minor
plot3(xnod(:,X), xnod(:,Y), xnod(:,Z), 'r*');
text(xnod(:,X), xnod(:,Y), xnod(:,Z), num2str((1:nno)'));
title('Numeracion de la estructura');
xlabel('Eje X [mm]')
ylabel('Eje Y [mm]')
zlabel('Eje Z [mm]')
box off

%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
c = [1:36]';              % los nodos 1, 2, 3, 4, 5, 6 estan empotrados
d = setdiff(1:ngdl, c);

%% ensamblo la matriz de rigidez global
K   = sparse(ngdl,ngdl);   % separo memoria
idx = cell(nbar,1);
T   = cell(nbar,1);
for e = 1:nbar  % para cada barra
   idx{e} = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) ]; % saco los 12 gdls de la barra
   
   V1 = xnod(LaG(e,1), :);
   V2 = xnod(LaG(e,2), :);
   V2_V1 = V2 - V1;
   L = norm(V2_V1);
%  V3_V1 = V3(e,:) - V1;
   e1g = (V2_V1)/L;
   e3g = cross(V2_V1,V3_V1(e,:)); e3g = e3g/norm(e3g);
   e2g = cross(e3g, e1g);
   
   TT = [ e1g; e2g; e3g ];      
   T{e} = blkdiag(TT, TT, TT, TT);

   AE = A(e)*E;   EIy = E*Iy(e);   EIz = E*Iz(e);   GJ = G*J(e); 
   L2 = L^2; L3 = L^3;
      
   Kloc = ...   
        [ AE/L,            0,            0,     0,           0,           0, -AE/L,            0,            0,     0,           0,           0
             0,  (12*EIz)/L3,            0,     0,           0,  (6*EIz)/L2,     0, -(12*EIz)/L3,            0,     0,           0,  (6*EIz)/L2
             0,            0,  (12*EIy)/L3,     0,  (6*EIy)/L2,           0,     0,            0, -(12*EIy)/L3,     0,  (6*EIy)/L2,           0
             0,            0,            0,  GJ/L,           0,           0,     0,            0,            0, -GJ/L,           0,           0
             0,            0,   (6*EIy)/L2,     0,   (4*EIy)/L,           0,     0,            0,  -(6*EIy)/L2,     0,   (2*EIy)/L,           0
             0,   (6*EIz)/L2,            0,     0,           0,   (4*EIz)/L,     0,  -(6*EIz)/L2,            0,     0,           0,   (2*EIz)/L
         -AE/L,            0,            0,     0,           0,           0,  AE/L,            0,            0,     0,           0,           0
             0, -(12*EIz)/L3,            0,     0,           0, -(6*EIz)/L2,     0,  (12*EIz)/L3,            0,     0,           0, -(6*EIz)/L2
             0,            0, -(12*EIy)/L3,     0, -(6*EIy)/L2,           0,     0,            0,  (12*EIy)/L3,     0, -(6*EIy)/L2,           0
             0,            0,            0, -GJ/L,           0,           0,     0,            0,            0,  GJ/L,           0,           0
             0,            0,   (6*EIy)/L2,     0,   (2*EIy)/L,           0,     0,            0,  -(6*EIy)/L2,     0,   (4*EIy)/L,           0
             0,   (6*EIz)/L2,            0,     0,           0,   (2*EIz)/L,     0,  -(6*EIz)/L2,            0,     0,           0,   (4*EIz)/L ];
    
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + T{e}'*Kloc*T{e}; % sumo a K global
end

%%
% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd   |   | Kcc Kcd || ac |   | fd |  Recuerde que qc = 0
%|      | = |         ||    | - |    |
%| qc=0 |   | Kdc Kdd || ad |   | fc | 

% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = zeros(length(c),1); % desplazamientos conocidos en contorno

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);        % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a  = zeros(ngdl,1);   a(c) = ac;  a(d) = ad; % desplazamientos
q  = zeros(ngdl,1);   q(c) = qd;             % fuerzas nodales equivalentes

%{

%% calculo las cargas axiales (N) en cada barra
axial = zeros(nbar,1);
for e = 1:nbar % para cada barra
   axial(e) = Taxial{e}*a(idx{e});
end;

%% imprimo los resultados
format short g
disp('Desplazamientos nodales                                                ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
vect_mov = reshape(a,3,nno)'; % vector de movimientos
for i = 1:nno
   fprintf('Nodo %3d: u = %12.4g mm, v = %12.4g mm, w = %12.4g mm\n', ...
      i, vect_mov(i,X), vect_mov(i,Y), vect_mov(i,Z));
end;

disp(' ');
disp('Fuerzas nodales de equilibrio (solo se imprimen los diferentes de cero)  ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
q = reshape(q,3,nno)';
for i = 1:nno   
   if ~isequal(q(i,:),[0 0 0])
      fprintf('Nodo %3d qx = %12.4g kN, qy = %12.4g kN, qz = %12.4g kN\n', ...
         i, q(i,X), q(i,Y), q(i,Z));
   end;
end;

disp(' ');
disp('Fuerzas axiales                    ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
for e = 1:nbar
   fprintf('Barra %3d: Axial = %12.4g kN, Esfuerzo = %12.4g kN/mm2\n', e, axial(e), axial(e)/A);
end;
%}

%% Dibujar la estructura y su deformada
esc.def    = 50;          % escalamiento de la deformada
esc.faxial = 0.2;         % escalamiento del diagrama de axiales
esc.V      = 0.3;         % escalamiento del diagrama de cortantes
esc.M      = 0.3;         % escalamiento del diagrama de momentos flectores
esc.T      = 0.3;         % escalamiento del diagrama de momentos torsores

%{
figure
hold on
view(3)
vect_mov = reshape(a,6,nno)'; % vector de movimientos
xdef = xnod + esc.def*vect_mov(:,[X Y Z]);
for e = 1:nbar
   plot3(xnod(LaG(e,1:2),X), xnod(LaG(e,1:2),Y), xnod(LaG(e,1:2),Z), 'b-');
   plot3(xdef(LaG(e,1:2),X), xdef(LaG(e,1:2),Y), xdef(LaG(e,1:2),Z), 'r-', 'LineWidth', 2);
end
xlabel('x, mm')
ylabel('y, mm')
zlabel('z, mm')
axis equal
grid minor
%}
%% bye, bye!
%return;



%% Dibujar la estructura y su deformada
figure(2); hold on; title('Deformada exagerada');       xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); axis equal
% figure(3); hold on; title('Fuerza axial [kN]');         xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); axis equal
% figure(4); hold on; title('Fuerza cortante Vy [kN]');   xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); axis equal
% figure(5); hold on; title('Fuerza cortante Vz [kN]');   xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); axis equal
% figure(6); hold on; title('Momento flector My [kN-m]'); xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); axis equal
% figure(7); hold on; title('Momento flector Mz [kN-m]'); xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); axis equal
% figure(8); hold on; title('Momento torsor T [kN-m]');   xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); axis equal

for e = 1:nbar
   V1 = xnod(LaG(e,1), :);  V2 = xnod(LaG(e,2), :);  
   material.A  = A(e);      material.E  = E;
   material.Iy = Iy(e);     material.Iz = Iz(e);
      
   carga.qxloc = @(x) 0; %qxloc{e};  % axial
   carga.qyloc = @(x) 0; %qyloc{e};  % en dir y
   carga.qzloc = @(x) 0; %qzloc{e};  % en dir z
   carga.qtloc = @(x) 0; %qtloc{e};  % torsion
  
   qe_loc = T{e}*q(idx{e});   ae_loc = T{e}*a(idx{e});
   
   matT = T{e}(1:3,1:3)';
   
   c1_dibujar_barra_deformada_portico_3d(material, V1, V2, carga, ...
                                                qe_loc, ae_loc, matT, esc);
end
