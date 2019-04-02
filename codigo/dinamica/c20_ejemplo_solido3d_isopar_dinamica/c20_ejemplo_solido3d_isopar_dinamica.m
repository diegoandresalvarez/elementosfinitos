clear, clc, close all; % borro la memoria, la pantalla y cierro los graficos
opengl software;       % para que funcione en Linux (o sino hace un segment fault 
                       % cuando hace los graficos)

%% defino las constantes
X = 1; Y = 2; Z = 3; consistent = 1; lumped = 2;
rhoe  = 7.860;        % densidad del material (ton/m^3)
Ee   = 200e6;         % modulo de elasticidad del sólido (kPa) = 200GPa
nue  = 0.33;          % coeficiente de Poisson
g    = 9.81;          % aceleracion de la gravedad (m/s^2)
be = [0; 0; -rhoe*g]; % vector de fuerzas másicas del elemento

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
ngdl = 3*nno;        % numero de grados de libertad
gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad
% gdl fila = nodo
% [gdl col1, gdl col2, gdl col3] = gdl en direccion [x, y z]

%% Se dibuja la malla de elementos finitos
figure;
tetramesh(LaG, xnod, 'FaceColor', 'cyan');
title('Malla de elementos finitos','FontSize',26);
daspect([1 1 1]); % similar a axis equal, pero en 3D
grid on;          % mostrar la rejilla
maximize;

%% Relacion de cargas puntuales
f = zeros(ngdl,1); % vector de fuerzas nodales equivalentes global
%f(gdl(4,Y)) = -1200;  % N - carga puntual en el nodo 4 dir Y
%f(gdl(5,X)) =   400;  % N - carga puntual en el nodo 5 dir X

%% Funciones de forma serendipitas del elemento tetraedrico de 4 nodos:
Nforma = @(alpha,beta,gamma) [ ...
1 - alpha - beta - gamma         % N1
alpha                            % N2
beta                             % N3
gamma                        ];  % N4

%% Derivadas de N con respecto a alpha
dN_dalpha = @(alpha,beta,gamma) [ ...
-1                               % dN1_dalpha
 1                               % dN2_dalpha
 0                               % dN3_dalpha
 0                           ];  % dN4_dalpha

%% Derivadas de N con respecto a beta
dN_dbeta = @(alpha,beta,gamma) [ ...
-1                               % dN1_dbeta
 0                               % dN2_dbeta
 1                               % dN3_dbeta
 0                           ];  % dN4_dbeta

%% Derivadas de N con respecto a gamma
dN_dgamma = @(alpha,beta,gamma) [ ...
-1                               % dN1_dgamma
 0                               % dN2_dgamma
 0                               % dN3_dgamma
 1                           ];  % dN4_dgamma

%% Parametros de la cuadratura de Gauss-Legendre
np = 4; % 1, 4 o 5 (numero de puntos de la cuadratura)
[L2L3L4_gl, w_gl] = tetra_quad(np);

%% ensamblo la matriz de rigidez K y matriz de masa M global
M = sparse(ngdl,ngdl); % matriz de masa global como RALA (sparse)
K = sparse(ngdl,ngdl); % matriz de rigidez global como RALA (sparse)
N = cell(nef,1);       % contenedor para las matrices de forma
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
   % Calculo las matrices de masa, rigidez y el vector de fuerzas nodales
   % equivalentes del elemento
   Me = zeros(12);
   Ke = zeros(12);
   fe = zeros(12,1);
   MATdet_Je = zeros(np,1); % en este vector se almacenarán los Jacobianos
   
   for p = 1:np   
      alpha_gl = L2L3L4_gl(p,1);
      beta_gl  = L2L3L4_gl(p,2);
      gamma_gl = L2L3L4_gl(p,3);
   
      % Se evaluan las funciones de forma en los puntos de integracion
      % de Gauss-Legendre
      NNforma = Nforma(alpha_gl, beta_gl, gamma_gl);
   
      % Se evaluan las derivadas de las funciones de forma en los puntos
      % de integracion de Gauss-Legendre
      ddN_dalpha = dN_dalpha(alpha_gl, beta_gl, gamma_gl); xe = xnod(LaG(e,:),X);
      ddN_dbeta  = dN_dbeta (alpha_gl, beta_gl, gamma_gl); ye = xnod(LaG(e,:),Y);
      ddN_dgamma = dN_dgamma(alpha_gl, beta_gl, gamma_gl); ze = xnod(LaG(e,:),Z);
         
      dx_dalpha = sum(ddN_dalpha.* xe);  dy_dalpha = sum(ddN_dalpha.* ye);   dz_dalpha = sum(ddN_dalpha.* ze);
      dx_dbeta  = sum(ddN_dbeta .* xe);  dy_dbeta  = sum(ddN_dbeta .* ye);   dz_dbeta  = sum(ddN_dbeta .* ze);
      dx_dgamma = sum(ddN_dgamma.* xe);  dy_dgamma = sum(ddN_dgamma.* ye);   dz_dgamma = sum(ddN_dgamma.* ze);      
     
      % Se ensambla la matriz Jacobiana del elemento
      Je = [ dx_dalpha dy_dalpha dz_dalpha
             dx_dbeta  dy_dbeta  dz_dbeta
             dx_dgamma dy_dgamma dz_dgamma ];
            
      % Se calcula el determinante del Jacobiano
      det_Je = det(Je);
      MATdet_Je(p) = det_Je;
      
      N{e} = zeros(3,3*4);
      B{e} = zeros(6,3*4);
      for i = 1:4
         % Se ensambla la matriz de funciones de forma N
         N{e}(:,[3*i-2 3*i-1 3*i]) = [ NNforma(i)  0           0         
                                       0           NNforma(i)  0
                                       0           0           NNforma(i) ];
      
         % Se ensambla la matriz de deformacion del elemento B
         tmp = Je\[ddN_dalpha(i); ddN_dbeta(i); ddN_dgamma(i)];         
         dNi_dx = tmp(1);
         dNi_dy = tmp(2);
         dNi_dz = tmp(3);
         
         % aquí se ensambla y asigna la matriz B_i
         B{e}(:,[3*i-2 3*i-1 3*i]) = [ dNi_dx 0      0  
                                       0      dNi_dy 0
                                       0      0      dNi_dz
                                       dNi_dy dNi_dx 0
                                       dNi_dz 0      dNi_dx
                                       0      dNi_dz dNi_dy ];
      end;
      
      % se arma la matriz de masa del elemento e
      switch mass_matrix
         case consistent
            Me = Me + rhoe*N{e}'*N{e}*det_Je*w_gl(p);
         case lumped
            Me = NaN; %no se ha definido aun
         otherwise
            error('La matriz de masa debe ser lumped o consistent');
      end;
      
      % se arma la matriz de rigidez del elemento e
      Ke = Ke + B{e}'*De*B{e}*det_Je*w_gl(p);

      % vector de fuerzas nodales equivalentes
      fe = fe + N{e}'*be*det_Je*w_gl(p);
   end
   
   if any(MATdet_Je <= 0)
      error('Existen elementos con det_Je negativo en el elemento %d.\n', e);
   end;
   
   % saco los 12 gdls del tetraedro
   idx = [gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) gdl(LaG(e,4),:)]; 
   
   M(idx,idx) = M(idx,idx) + Me;   
   K(idx,idx) = K(idx,idx) + Ke;
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
ad = Kdd\(fc - Kdc*ac);    % calculo desplazamientos desconocidos
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
escala = 5e3; % factor de escalamiento de la deformada
xdef = xnod + escala*delta; % posicion de la deformada
figure
tetramesh(LaG, xdef, 'FaceColor', 'blue');
title(sprintf('Deformada escalada %d veces',escala));
daspect([1 1 1]); % similar a axis equal, pero en 3D
grid on;          % mostrar la rejilla
maximize;

%% **** **** **** ANALISIS DINAMICO **** **** ****
% estimamos las frecuencias naturales y los modos de vibracion del solido

% Calcula los modos de vibracion y las frecuencias naturales de una 
% estructura con matriz de rigidez K y matriz de masa M
n = 20; % numero de modos a tener en cuenta en el calculo

OPT.maxit = 10000;
OPT.tol   = 1e-32;
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
disp('Frecuencias de vibracion');
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
escala = 0.2;   % ayuda a la visualización
figure         % crea un lienzo
hold on;       % permite que se sobreescriba sobre dicho lienzo
maximize       % maximize la ventana
axis([-1.5 1.5 -1.5 1.5 0 1.5]); % ejes 
grid on        % cuadricula
%tetramesh(LaG, xnod, 'FaceColor', 'cyan','EdgeAlpha', 0.1,'FaceAlpha',0.1);
daspect([1 1 1]); % similar a axis equal, pero en 3D
view(3)

ncuadros = 50; % numero de cuadros de la animacion
t = linspace(0,2*pi,ncuadros); % intervalos de tiempo en la animacion

for mv = 1:1   
   title(sprintf('Modo %3d   w = %12.4f rad,   f = %12.4f Hz,   T = %12.4f seg \n', ...
      mv, Omega(mv), 1/T(mv), T(mv)), 'FontSize', 15);
   
   def = escala*reshape(PHI(:,mv)',3,nno)';

% w = waitbar(0,'My Progress Bar');   
   for i = 2:ncuadros
      fprintf('mv = %d, i = %d\n',mv,i);
%     w = waitbar((i-1)/(ncuadros-1),w,['Iteracion: ',num2str(i-1)]);     % update the wait bar each iteration      
      xdef = xnod + def*sin(t(i)); % modo de vibracion
%      deftot = sqrt(def(:,1).^2 + def(:,2).^2 + def(:,3).^2);
      trep = TriRep(LaG, xdef);
      [tri,xf] = freeBoundary(trep);
      h = trisurf(tri, xf(:,1),xf(:,2),xf(:,3), 'FaceColor', 'red'); % 'EdgeAlpha', 0.1,'FaceAlpha',0.1);
      
      drawnow; % grafique todo lo que tenga pendiente antes de continuar
      print('-dpng', '-zbuffer', '-r200', sprintf('xsol%02d_%02d.png', mv,i));
      delete(h); % borra las deformada contenida en h
   end

   % converti los .png a .gif con el comando de linux:
   unix(sprintf('convert -antialias -loop 10000 -delay 5 -scale 80%% xsol* c20_ejemplo_solido3d_mv_%02d.gif',mv));
   unix(sprintf('rm xsol%02d_*.png',mv));
   %web('c20_ejemplo_solido3d_modos_vibracion.gif');
%   close(w)
end
close
