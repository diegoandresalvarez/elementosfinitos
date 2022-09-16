clear, clc, close all  % borro la memoria, la pantalla y cierro los graficos
%opengl software;       % para que funcione en Linux (o sino hace un segment fault 
                       % cuando hace los graficos)

%% defino las constantes
X = 1; Y = 2; Z = 3; 
consistent = 1; lumped = 2;
rhoe  = 7860;       % kg/m^3 densidad del material
Ee   = 200e9;       % modulo de elasticidad del sólido (Pa) = 200GPa
nue  = 0.33;        % coeficiente de Poisson
g    = 9.81;        % aceleracion de la gravedad (m/s^2)

%% defino las variables
mass_matrix = consistent;         % matriz de masa lumped o consistent ?

%% Cargo en la memoria las variables xnod, LaG y nod_ac
solido3d;
% este archivo define 
% xnod   = definicion de las coordenadas nodales (x,y,z) en m
% nod_ac = nodos con desplazamiento conocido igual a cero (apoyos)

% LaG    = definicion de los tetraedros (i,j,k,l)
% LaG: local a global: matriz que relaciona gdls locales y globales
%     fila = barra
%     col1 = gdl global asociado a gdl local 1,
%     col2 = gdl global asociado a gdl local 2,
%     col3 = gdl global asociado a gdl local 3
%     (se lee el tetraedro en la fila x tiene vertices 1, 2 y 3


%% Numero de nodos, elementos finitos y grados de libertad
nno  = size(xnod,1); % numero de nodos
nef  = size(LaG,1);  % numero de elementos finitos (numero de tetraedros)
ngdl = 3*nno; % numero de grados de libertad
gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad
% gdl fila = nodo
% gdl col1 = gdl en direccion x
% gdl col2 = gdl en direccion y
% gdl col3 = gdl en direccion z                                


%% Se dibuja la malla de elementos finitos
figure;
tetramesh(LaG, xnod, 'FaceColor', 'cyan');
title('Malla de elementos finitos');
daspect([1 1 1]); % similar a axis equal, pero en 3D
grid on;          % mostrar la rejilla
maximize;

%% Relacion de cargas puntuales
f = zeros(ngdl,1); % vector de fuerzas nodales equivalentes global
%f(gdl(4,Y)) = -1200;  % N - carga puntual en el nodo 4 dir Y
%f(gdl(5,X)) =   400;  % N - carga puntual en el nodo 5 dir X

%% ensamblo la matriz de rigidez K y matriz de masa M global
M = sparse(ngdl,ngdl); % matriz de masa global como RALA (sparse)
K = sparse(ngdl,ngdl); % matriz de rigidez global como RALA (sparse)
B = cell(nef,1);       % contenedor para las matrices de deformación

% matriz constitutiva del elemento 3D
d1 = (1-nue)/(1-2*nue);
d2 = nue/(1-2*nue);
De = Ee/(1+nue)*[ d1 d2 d2 0   0   0
                  d2 d1 d2 0   0   0
                  d2 d2 d1 0   0   0
                  0  0  0  1/2 0   0
                  0  0  0  0   1/2 0
                  0  0  0  0   0   1/2 ];

for e = 1:nef    % ciclo sobre todos los elementos finitos  
   % Determino las coordenadas del los nodos y el volumen del EF
   x1 = xnod(LaG(e,1),X); y1 = xnod(LaG(e,1),Y); z1 = xnod(LaG(e,1),Z);
   x2 = xnod(LaG(e,2),X); y2 = xnod(LaG(e,2),Y); z2 = xnod(LaG(e,2),Z);
   x3 = xnod(LaG(e,3),X); y3 = xnod(LaG(e,3),Y); z3 = xnod(LaG(e,3),Z);
   x4 = xnod(LaG(e,4),X); y4 = xnod(LaG(e,4),Y); z4 = xnod(LaG(e,4),Z);
   
   Ve = (1/6)*det([ 1 x1 y1 z1      %Volumen del EF e
                    1 x2 y2 z2
                    1 x3 y3 z3
                    1 x4 y4 z4]);               
   if Ve <= 0
      error(sprinft('Revise las coordenadas locales del EF %d.\n', e));
   end
   
   % Calculo de la matriz de masa del elemento e
   switch mass_matrix
      case consistent  % matriz de masa CONSISTENTE
         Me = (rhoe*Ve/20)*...
             [ 2 0 0 1 0 0 1 0 0 1 0 0
               0 2 0 0 1 0 0 1 0 0 1 0
               0 0 2 0 0 1 0 0 1 0 0 1
               1 0 0 2 0 0 1 0 0 1 0 0
               0 1 0 0 2 0 0 1 0 0 1 0
               0 0 1 0 0 2 0 0 1 0 0 1
               1 0 0 1 0 0 2 0 0 1 0 0
               0 1 0 0 1 0 0 2 0 0 1 0
               0 0 1 0 0 1 0 0 2 0 0 1
               1 0 0 1 0 0 1 0 0 2 0 0
               0 1 0 0 1 0 0 1 0 0 2 0
               0 0 1 0 0 1 0 0 1 0 0 2 ];
      case lumped % matriz de masa local CONCENTRADA         
         Me = (rhoe*Ve/4)*diag(ones(12,1));
      otherwise
         error('Tipo de matriz de masa invalido');
   end   
   
   % Cálculo de la matriz de deformaciones B.
   a1 =   det([x2 y2 z2; x3 y3 z3; x4 y4 z4]);
   a2 =  -det([x3 y3 z3; x4 y4 z4; x1 y1 z1]);
   a3 =   det([x4 y4 z4; x1 y1 z1; x2 y2 z2]);
   a4 =  -det([x1 y1 z1; x2 y2 z2; x3 y3 z3]);
   
   b1 = -det([ 1 y2 z2;  1 y3 z3;  1 y4 z4]); % OJO CON LOS SIGNOS
   b2 =  det([ 1 y3 z3;  1 y4 z4;  1 y1 z1]); % NO DAN LO MISMO QUE LA
   b3 = -det([ 1 y4 z4;  1 y1 z1;  1 y2 z2]); % LITERATURA
   b4 =  det([ 1 y1 z1;  1 y2 z2;  1 y3 z3]);
   
   c1 = -det([x2  1 z2; x3  1 z3; x4  1 z4]);
   c2 =  det([x3  1 z3; x4  1 z4; x1  1 z1]);
   c3 = -det([x4  1 z4; x1  1 z1; x2  1 z2]);
   c4 =  det([x1  1 z1; x2  1 z2; x3  1 z3]);
   
   d1 = -det([x2 y2  1; x3 y3  1; x4 y4  1]);
   d2 =  det([x3 y3  1; x4 y4  1; x1 y1  1]);
   d3 = -det([x4 y4  1; x1 y1  1; x2 y2  1]);
   d4 =  det([x1 y1  1; x2 y2  1; x3 y3  1]);
   
   B{e} = (1/(6*Ve))*[ b1  0  0   b2  0  0   b3  0  0   b4  0  0
                        0 c1  0    0 c2  0    0 c3  0    0 c4  0
                        0  0 d1    0  0 d2    0  0 d3    0  0 d4
                       c1 b1  0   c2 b2  0   c3 b3  0   c4 b4  0
                       d1  0 b1   d2  0 b2   d3  0 b3   d4  0 b4
                       0  d1 c1   0  d2 c2   0  d3 c3   0  d4 c4 ];

   % Calculo de la matriz de rigidez del elemento e
   Ke = B{e}'*De*B{e}*Ve;
              
   % Calculo del vector de fuerzas nodales equivalentes del elemento e
   % Fuerzas masicas (peso propio)
   fbe = [0; 0; -1; 
          0; 0; -1;
          0; 0; -1;
          0; 0; -1] * rhoe*g*Ve/4;
        
   fe = fbe; % vector de fuerzas nodales equivalentes
   
   % saco los 12 gdls del tetraedro
   idx = [gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) gdl(LaG(e,4),:)]; 
   
   % Ensamblo las contribuciones a las matrices globales   
   M(idx,idx) = M(idx,idx) + Me; % sumo a M global   
   K(idx,idx) = K(idx,idx) + Ke; % sumo a K global
   f(idx,:)   = f(idx,:)   + fe;
   % NOTA: para quitar el Warning de MATLAB que el indexing is slow mire:
   % http://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finite-element-matrices-in-matlab/
end;

if mass_matrix == lumped % se remueven los errores de redondeo que aparecen
   M = diag(diag(M));    % fuera de la diagonal. Mire la ayuda de funcion diag
end;

%% Muestro la configuración de las matrices M y K (ambas son ralas)
figure
subplot(1,2,1); spy(M); title('Matriz de masa M global');
subplot(1,2,2); spy(K); title('Matriz de rigidez K global');
maximize;

%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
c  = [gdl(nod_ac,X); gdl(nod_ac,Y); gdl(nod_ac,Z)]; d = setdiff(1:ngdl,c);
ac = zeros(length(c),1); % desplazamientos conocidos

%| Mcc Mcd ||ddac|   | Kcc Kcd || ac |   | fd |   |  qd  |
%|         ||ddac| + |         ||    | - |    | = |      |
%| Mdc Mdd ||ddac|   | Kdc Kdd || ad |   | fc |   | qc=0 |

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%% **** **** **** ANALISIS ESTATICO **** **** ****
% extraigo las submatrices y especifico las cantidades conocidas
Mcc = M(c,c); Mcd = M(c,d); Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Mdc = M(d,c); Mdd = M(d,d); Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% resuelvo el sistema de ecuaciones
%ad = Kdd\(fc - Kdc*ac);    % calculo desplazamientos desconocidos
AA = mean(mean(abs(Kdd)));
bb = mean(mean(abs(fc - Kdc*ac)));

ad = (Kdd/AA)\((fc - Kdc*ac)/bb);    % calculo desplazamientos desconocidos
ad = ad*bb/AA;

qd = Kcc*ac + Kcd*ad - fd; % calculo fuerzas de equilibrio desconocidas
a = zeros(ngdl,1);  a(c) = ac;  a(d) = ad; % desplazamientos
q = zeros(ngdl,1);  q(c) = qd; %q(d) = qc; % fuerzas nodales equivalentes

%% imprimo los resultados
disp('Desplazamientos nodales                      ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
a = reshape(a,3,nno)';
for i = 1:nno
   fprintf('Nodo %3d u = %12.4g m, v = %12.4g m, w = %12.4g m\n', i, a(i,X), a(i,Y), a(i,Z));
end;

disp(' ');
disp('Fuerzas nodales de equilibrio (solo se imprimen los diferentes de cero)');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
q = reshape(q,3,nno)';
for i = 1:nno   
   if ~isequal(q(i,:),[0 0 0])
      fprintf('Nodo %3d qx = %12.4g N, qy = %12.4g N, qz = %12.4g N\n', i, q(i,X), q(i,Y), q(i,Z));
   end;
end;

%% Dibujo la malla de elementos finitos y las deformaciones de esta
delta = reshape(a,3,nno)';
escala = 1000; % factor de escalamiento de la deformada
xdef = xnod + escala*delta; % posicion de la deformada
figure
%tetramesh(LaG, xnod, 'FaceColor', 'cyan');
tetramesh(LaG, xdef, 'FaceColor', 'blue');
%hold on 
title(sprintf('Deformada escalada %d veces',escala));
daspect([1 1 1]); % similar a axis equal, pero en 3D
grid on;          % mostrar la rejilla
maximize;

%% **** **** **** ANALISIS DINAMICO **** **** ****
% estimamos las frecuencias naturales y los modos de vibracion del solido

% Calcula los modos de vibracion y las frecuencias naturales de una 
% estructura con matriz de rigidez K y matriz de masa M
n = 15; % numero de modos a tener en cuenta en el calculo

OPT.maxit = 1000;
OPT.tol   = 1e-16;
% OPT.issym = 1;  % La matriz Kdd es simetrica
[Phi,Omega2,flag] = eigs(Kdd,Mdd,n,'SM',OPT);% autovalores y autovectores
if flag ~= 0
   warning('El calculo de valores y vectores propios no convergio');
end
%{
tic
[Phi,Omega2] = eig(full(Kdd),full(Mdd));% autovalores y autovectores
toc
%}

Omega = sqrt(diag(Omega2));  % frecuencias angulares
[Omega,I] = sort(Omega);     % ordena las frecuencias
T = 2*pi./Omega;             % seg - periodo de vibracion 
Phi = Phi(:,I); % ordena los modos segun el orden de las frecuencias
Phi = Phi./repmat(sqrt(diag(Phi'*Mdd*Phi)'),size(Phi,1),1); % normaliza los modos

PHI      = zeros(length(d),n);
PHI(d,:) = Phi(:,1:n);
PHI(c,:) = repmat(ac,1,n);

disp(' ');
disp('Frecuencias de vibracion                     ');
switch mass_matrix
   case consistent, disp('(para una matriz de masa consistente)');
   case lumped,     disp('(para una matriz de masa condensada)'); 
end
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
for i = 1:n
   fprintf('Modo %3d   w = %12.4f rad,   f = %12.4f Hz,   T = %12.4f seg \n', ...
      i, Omega(i), 1/T(i), T(i));
end;

%% Hacer animacion del modo de vibracion "mv"
escala = 10;   % ayuda a la visualización
figure         % crea un lienzo
hold on;       % permite que se sobreescriba sobre dicho lienzo
maximize       % maximize la ventana
axis([-1.5 1.5 -1.5 1.5 0 1.5]); % ejes 
grid on        % cuadricula
tetramesh(LaG, xnod, 'FaceColor', 'cyan');
daspect([1 1 1]); % similar a axis equal, pero en 3D
view(3)

ncuadros = 50; % numero de cuadros de la animacion
t = linspace(0,2*pi,ncuadros); % intervalos de tiempo en la animacion

for mv = 1:1;   %1:Nummodos_a_graficar
   title(sprintf('Modo %3d   w = %12.4f rad,   f = %12.4f Hz,   T = %12.4f seg \n', ...
      mv, Omega(mv), 1/T(mv), T(mv)), 'FontSize', 15);
   
   def = escala*reshape(PHI(:,mv)',3,nno)';

   % w = waitbar(0,'My Progress Bar');   
   for i = 2:ncuadros
      % w = waitbar((i-1)/(ncuadros-1),w,['Iteracion: ',num2str(i-1)]);     % update the wait bar each iteration      
      xdef = xnod + def*sin(t(i)); % modo de vibracion
      h = tetramesh(LaG, xdef, 'FaceColor', 'red');
      drawnow; % grafique todo lo que tenga pendiente antes de continuar
      print('-dpng', '-zbuffer', '-r200', sprintf('xsol%02d.png', i));
      delete(h); % borra las deformada contenida en h
   end
   % close(w)
end
close

% converti los .png a .gif con el comando de linux:
unix('convert -antialias -loop 10000 -delay 5 -scale 80% xsol* c20_ejemplo_solido3d_modos_vibracion.gif');
unix('rm xsol*.png');
web('c20_ejemplo_solido3d_modos_vibracion.gif');
