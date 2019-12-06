clear, clc, close all

%% constantes
X = 1; Y = 2; Z = 3;

%% Unidades en kN y mm

% Cargar xnod y LaG
% xnod: coordenadas de cada nodo [x, y, z]
% LaG: local a global: matriz que relaciona nodos locales y globales
% fila = barra
% col1 = nodo global asociado a nodo local 1
% col2 = nodo global asociado a nodo local 2   
% (se lee la barra x va del nodo i al nodo j)
load domo;

A = 100;    % area barra (mm^2)
E = 205.8;  % modulo de elasticidad (kN/mm^2)

nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
nbar = size(LaG,1);  % numero de EFs (numero de filas de LaG)
ngdl = 3*nno;        % numero de grados de libertad (tres por nodo)

%% gdl: grados de libertad
% fila = nodo
% col1 = gdl en direccion x
% col2 = gdl en direccion y
% col3 = gdl en direccion z
gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad

%% cargas aplicadas
nodos_cargados = [1:10 16:21 26:31 36:41 46:51 56:58];
dofs_cargados  = 3*nodos_cargados;  % cargas en la direccion z

f = zeros(ngdl, 1);
f(dofs_cargados) = -50/3; %kN


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

%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
c = [31:45 64:75 94:105 124:135 154:165 175:183]';
d = setdiff(1:ngdl, c);

%% ensamblo la matriz de rigidez global
K = zeros(ngdl);   % separo memoria
Taxial = cell(nbar,1);
for e = 1:nbar  % para cada barra
   idx = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) ]; % saco los 6 gdls de la barra
   
   %{
     dxdydz = xnod(LaG(e,2), :) - xnod(LaG(e,1), :);
     L = norm(dxdydz); 
     l = dxdydz(1)/L;
     m = dxdydz(2)/L;
     n = dxdydz(3)/L;
   %}
   
   x1 = xnod(LaG(e,1), X);  x2 = xnod(LaG(e,2), X);
   y1 = xnod(LaG(e,1), Y);  y2 = xnod(LaG(e,2), Y);
   z1 = xnod(LaG(e,1), Z);  z2 = xnod(LaG(e,2), Z);
   
   L = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2);
   l = (x2-x1)/L;
   m = (y2-y1)/L;
   n = (z2-z1)/L;

   TT = [ l*l l*m l*n
          m*l m*m m*n
          n*l n*m n*n ]; % = [l m n].'*[l m n]

       
	kbar = A*E/L;
   
   Taxial{e} = kbar*[ -l -m -n l m n ];
      
   %{   
   syms l m n

   T = [l m n 0 0 0
        0 0 0 l m n];

   Kloc = [ 1 -1
           -1  1 ];

   K = T.'*Kloc*T
   %}
      
   K(idx,idx) = K(idx,idx) + kbar*[TT -TT; -TT TT]; % sumo a K global
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

%% calculo las cargas axiales (N) en cada barra
axial = zeros(nbar,1);
for e = 1:nbar % para cada barra
   idx = [ gdl(LaG(e,1),:) gdl(LaG(e,2),:) ]; % saco los 6 gdls de la barra
   axial(e) = Taxial{e}*a(idx);
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

%% Dibujar el domo, su deformada y las cargas

figure
hold on
view(3)
xdef = xnod + 50*vect_mov;
for e = 1:nbar
   plot3(xnod(LaG(e,1:2),X), xnod(LaG(e,1:2),Y), xnod(LaG(e,1:2),Z), 'b-');
   plot3(xdef(LaG(e,1:2),X), xdef(LaG(e,1:2),Y), xdef(LaG(e,1:2),Z), 'r-', 'LineWidth', 2);
end
xlabel('x, mm')
ylabel('y, mm')
zlabel('z, mm')
axis equal
grid minor
cargas_punt = 100*reshape(f,3,nno)'; % vector de movimientos
for i = 1:nno
   quiver3(xnod(i,X), xnod(i,Y), xnod(i,Z), ...
      cargas_punt(i,X), cargas_punt(i,Y), cargas_punt(i,Z), ...
      'k-', 'LineWidth',2, 'MaxHeadSize',0.5);
end

%%
return;

