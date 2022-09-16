clear, clc, close all % borro la memoria, la pantalla y cierro los graficos

%% defino las constantes
X = 1; Y = 2; consistent = 1; lumped = 2;

%% defino las variables
area  = repmat(2.5e-3, 1, 9); % m^2 area de las barra
rhoe  = 7860;                 % kg/m^3 densidad del material
E     = 200e9;                % Pa modulo de elasticidad del material
mass_matrix = lumped;         % matriz de masa lumped o consistent ?

%% Posicion (x,y) de cada nodo en cm
xnod = [ 0     0
         4     0
         4     3
         8     0
         8     3
        12     0 ];
nno = size(xnod,1); % numero de nodos
      
%% LaG: local a global: matriz que relaciona gdls locales y globales
LaG = [1 2   % fila = barra
       1 3   % col1 = gdl global asociado a gdl local 1,
       2 3   % col2 = gdl global asociado a gdl local 2,
       2 4 
       3 4
       3 5
       4 5
       4 6
       5 6 ]; % (se lee la barra en la fila x va del nodo i al nodo j)    
nef = size(LaG,1); % numero de elementos finitos (numero de barras)

%% gdl: grados de libertad
gdl   = [(1:2:2*nno)' (2:2:2*nno)'];  % fila = nodo
                                  % col1 = gdl en direccion x
                                  % col2 = gdl en direccion y
                                 
ngdl = 2*nno; % numero de grados de libertad

%% Se dibuja la malla de elementos finitos
figure;
hold on;
cgx = zeros(1,nef); cgy = zeros(1,nef); % almacena el centro de gravedad
for e = 1:nef
   line(xnod(LaG(e,[1 2]),X), xnod(LaG(e,[1 2]),Y));
   
   % Calculo la posicion del centro de gravedad del triangulo
   cgx(e) = (xnod(LaG(e,1),X) + xnod(LaG(e,2),X))/2;
   cgy(e) = (xnod(LaG(e,1),Y) + xnod(LaG(e,2),Y))/2;
   h = text(cgx(e), cgy(e), num2str(e)); set(h,'Color', [1 0 0]);
end
plot(xnod(:,X), xnod(:,Y), 'r*');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
axis equal tight
title('Malla de elementos finitos');

%% Relacion de cargas puntuales
f = zeros(ngdl,1); % vector de fuerzas nodales equivalentes global
f(gdl(4,Y)) = -1200;  % N - carga puntual en el nodo 4 dir Y
f(gdl(5,X)) =   400;  % N - carga puntual en el nodo 5 dir X

%% Se calcula la longitud y el angulo de inclinacion de cada barra
long  = zeros(1,nef);
theta = zeros(1,nef);
for e = 1:nef
   dx = xnod(LaG(e,2),X) - xnod(LaG(e,1),X);
   dy = xnod(LaG(e,2),Y) - xnod(LaG(e,1),Y);
   long(e)  = sqrt(dx^2 + dy^2); % longitud de cada barra
   theta(e) = atan2(dy,dx);      % angulo de inclinacion de cada barra
end

k = E.*area./long;  % rigidez de cada barra

%% ensamblo la matriz de rigidez global
T = cell(nef,1); % separo la memoria
M = zeros(ngdl);
K = zeros(ngdl);
for e = 1:nef    % para cada barra
   idx = [gdl(LaG(e,1),:) gdl(LaG(e,2),:)]; % saco los 4 gdls de la barra
   c = cos(theta(e)); s = sin(theta(e));  % sin y cos de la inclinacion
   
   T{e} = [ c  s  0  0     % matriz de transformacion de coordenadas
           -s  c  0  0     % para la barra e
            0  0  c  s
            0  0 -s  c ];
         
   switch mass_matrix
      case consistent
         Mloc = (rhoe*area(e)*long(e)/6)*...
            [ 2  0  1  0   % matriz de masa local CONSISTENTE
              0  2  0  1   % expresada en el sistema de coordenadas
              1  0  2  0   % locales para la barra e
              0  1  0  2 ];
      case lumped
         Mloc = (rhoe*area(e)*long(e)/2)*...
            [ 1  0  0  0   % matriz de masa local CONCENTRADA
              0  1  0  0   % expresada en el sistema de coordenadas
              0  0  1  0   % locales para la barra e
              0  0  0  1 ];
      otherwise
         error('Tipo de matriz de masa invalido');
   end

   Kloc = k(e)*[ 1  0 -1  0   % matriz de rigidez local expresada en el 
                 0  0  0  0   % sistema de coordenadas locales para la                  
                -1  0  1  0   % barra e
                 0  0  0  0 ];
              
   % Calculo del vector de fuerzas nodales equivalentes del elemento e
   % Fuerzas masicas (peso propio): observe que ya est� dada en coordenadas
   % globales
   fbe = [0; -1; 0; -1] * rhoe*area(e)*long(e)/2;
        
   fe = fbe; % vector de fuerzas nodales equivalentes
   
   % Ensamblo las contribuciones a las matrices globales   
   M(idx,idx) = M(idx,idx) + T{e}'*Mloc*T{e}; % sumo a M global
   K(idx,idx) = K(idx,idx) + T{e}'*Kloc*T{e}; % sumo a K global   
   f(idx,:)   = f(idx,:)   + fe;   
end;

if mass_matrix == lumped % se remueven los errores de redondeo que aparecen
   M = diag(diag(M));    % fuera de la diagonal. Mire la ayuda de funcion diag
end;

%% matriz masa y de rigidez global
figure
subplot(1,2,1); spy(M); title('Matriz de masa M global');
subplot(1,2,2); spy(K); title('Matriz de rigidez K global');

%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
c  = [gdl(1,X) gdl(1,Y) gdl(6,Y)]; d = setdiff(1:ngdl,c);
ac = [0;       0;       0       ]; % desplaz. conocidos para gdl(1,X) gdl(1,Y) gdl(2,Y)

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

% calculo las cargas axiales (N) en cada barra
N = zeros(nef,1);
for e = 1:nef % para cada barra
   idx = [gdl(LaG(e,1),:) gdl(LaG(e,2),:)]; % saco los 4 gdls de la barra e
   N(e) = k(e)*[-1 0 1 0]*T{e}*a(idx);
end;

%% imprimo los resultados
disp('Desplazamientos nodales                      ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
a = reshape(a,2,nno)';
for i = 1:nno
   fprintf('Nodo %3d  u = %12.4g m,   v = %12.4g m \n', i, a(i,X), a(i,Y));
end;

disp(' ');
disp('Fuerzas nodales de equilibrio                ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
q = reshape(q,2,nno)';
for i = 1:nno
   fprintf('Nodo %3d qx = %12.4g N,  qy = %12.4g N \n', i, q(i,X), q(i,Y));
end;

disp(' ');
disp('Fuerzas axiales en la barra                  ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
for i = 1:nef
   fprintf('Barra %3d  N = %9.4g N \n', i, N(i));
end;

disp(' ');
disp('Esfuerzos en la barra                  ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
for i = 1:nef
   fprintf('Barra %3d  sigma = %9.4g Pa \n', i, N(i)/area(i));
end;


%% **** **** **** ANALISIS DINAMICO **** **** ****
% estimamos las frecuencias naturales y los modos de vibracion de la cercha

% Calcula los modos de vibracion y las frecuencias naturales de una 
% estructura con matriz de rigidez K y matriz de masa M
[Phi,Omega2] = eig(Kdd,Mdd); % autovectores y autovalores 
Omega = sqrt(diag(Omega2));  % frecuencias angulares
[Omega,I] = sort(Omega);     % ordena las frecuencias
T = 2*pi./Omega;             % seg - periodo de vibracion 
Phi = Phi(:,I); % ordena los modos segun el orden de las frecuencias
Phi = Phi./repmat(sqrt(diag(Phi'*Mdd*Phi)'),size(Phi,1),1); % normaliza los modos

n = length(d);        % n�mero de gdl a tener en cuenta en el c�lculo
PHI = zeros(ngdl,n);
PHI(d,:) = Phi;
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
axis([-1 13 -1 4]); % ejes 
grid on        % cuadricula
for e = 1:nef  % grafico cada barra de la cercha en color azul
   line(xnod(LaG(e,:),X), xnod(LaG(e,:),Y),'Color','b');
end

ncuadros = 50; % numero de cuadros de la animacion
t = linspace(0,2*pi,ncuadros); % intervalos de tiempo en la animacion
h = zeros(nef,1); % h contendra un indicador de las lineas del modo de vibracion

for mv = 1:n;   
   title(sprintf('Modo %3d   w = %12.4f rad,   f = %12.4f Hz,   T = %12.4f seg \n', ...
      mv, Omega(mv), 1/T(mv), T(mv)), 'FontSize', 20);
  
   def = escala*reshape(PHI(:,mv)',2,nno)';
  
   tic
   while toc < 3 % animar por 3 segundos
      for i = 2:ncuadros
         xdef = xnod + def*sin(t(i)); % modo de vibracion
         for e = 1:nef
            h(e) = line(xdef(LaG(e,:),X), xdef(LaG(e,:),Y), ...
               'Color','r', 'LineWidth', 5);
         end
         drawnow; % grafique todo lo que tenga pendiente antes de continuar
%         print('-dpng', '-zbuffer', '-r100', sprintf('cercha%02d.png', i));
         delete(h); % borra las lineas rojas contenidas en h
      end
   end
end
close

% converti los .png a .gif con el comando de linux:
% convert -antialias -loop 10000 -delay 5 -scale 80% cercha* c20_ejemplo_cercha_modos_vibracion.gif
