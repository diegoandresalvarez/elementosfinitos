% Ejemplo 11.23 Uribe Escamilla

%% Unidades en toneladas y metros

% se definen algunas constantes que hacen el código más legible
NL1 = 1; NL2 = 2;
X = 1;   Y = 2;

%% defino las variables
Aviga = 0.30*0.35;       Acol  = 0.30*0.30;       % m^2    area
Iviga = 0.30*0.35^3/12;  Icol  = 0.30*0.30^3/12;  % m^4    inercia_y

% barra   1             2          3
A     = [ Aviga         Acol       Acol          ]; % areas
I     = [ Iviga         Icol       Icol          ]; % inercias_y
long  = [ hypot(4,2)    5          hypot(2,6)    ]; % long barra (m)
theta = [ atan2(2,4)    atan2(4,3) atan2(-6,2)   ]*180/pi; % angulo inclinacion (grados)

% LaG: local a global: matriz que relaciona nodos locales y globales
% (se lee la barra x va del nodo i al nodo j)
LaG = [1 2    % fila = barra
       4 1    % col1 = nodo global asociado a nodo local 1
       2 3];  % col2 = nodo global asociado a nodo local 2

% gdl: grados de libertad
ngdl = 12; % numero de grados de libertad
gdl  = [ 4  5  6   % fila = nodo
         7  8  9   % col1 = gdl en direccion x
        10 11 12   % col2 = gdl en direccion y
         1  2  3]; % col3 = gdl en direccion angular antihoraria

E = 190*10000; %ton/m^2  modulo de elasticidad

nb = size(LaG,1); %numero de barras (numero de filas de LaG)
nn = size(gdl,1); %numero de nodos  (numero de filas de gdl)

%% fuerzas nodales equivalentes para las diferentes barras
% (en este ejemplo las fuerzas nodales equivalentes estas siendo 
% especificadas con respecto al sistema de coordenadas globales)
fe = cell(nb,1);
%        fxi    fyi    mi     fxj    fyj    mj
%        ton    ton    ton-m  ton    ton    ton-m
fe{1} = [0     -5.60  -3.733  0     -5.60   +3.733 ]'; % OJO con los signos
fe{2} = [0      0      0      0      0      0      ]'; % mirar pag 613
fe{3} = [0      0      0      0      0      0      ]';

%% separo la memoria
K   = zeros(ngdl);    % matriz de rigidez global
f   = zeros(ngdl,1);  % vector de fuerzas nodales equivalentes global
Ke  = cell(nb,1);     % matriz de rigidez local en coordenadas globales
T   = cell(nb,1);     % matriz de transformacion de coordenadas
idx = zeros(nb,6);    % almacena los 6 gdls de las barras

%% ensamblo la matriz de rigidez global (K) y vector de fuerzas global (f)
for e = 1:nb  % para cada barra
   % saco los 6 gdls de la barra e
   idx(e,:) = [gdl(LaG(e,NL1),:) gdl(LaG(e,NL2),:)];
   
   % matriz de transformacion de coordenadas para la barra e
   c = cosd(theta(e)); s = sind(theta(e));  % sin y cos de la inclinacion
   T{e} = [ c  s  0  0  0  0;        
           -s  c  0  0  0  0;        
            0  0  1  0  0  0;
            0  0  0  c  s  0;
            0  0  0 -s  c  0;
            0  0  0  0  0  1];
         
   % matriz de rigidez local expresada en el sistema de coordenadas locales
   % para la barra e
   AE = A(e)*E;    EI = E*I(e);    L=long(e); L2=long(e)^2; L3=long(e)^3;
   Kloc = [ AE/L   0         0        -AE/L    0          0       
            0     12*EI/L3   6*EI/L2   0     -12*EI/L3   6*EI/L2
            0      6*EI/L2   4*EI/L    0      -6*EI/L2   2*EI/L
           -AE/L   0         0         AE/L    0         0
            0    -12*EI/L3  -6*EI/L2   0      12*EI/L3  -6*EI/L2
            0      6*EI/L2   2*EI/L    0      -6*EI/L2   4*EI/L];

   % matriz de rigidez local en coordenadas globales
   Ke{e} = T{e}'*Kloc*T{e};            
   K(idx(e,:),idx(e,:)) = K(idx(e,:),idx(e,:)) + Ke{e}; % sumo Ke{e} a K global
   f(idx(e,:))          = f(idx(e,:))          + fe{e}; % sumo a f global
end;

%% localizo la carga puntual de 1.5 ton en el gdl 4
f(4) = f(4) + 1.5;

%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
c = [1 2 3 10 11 12];    
d = setdiff(1:ngdl, c); % d = [4 5 6 7 8 9];

%% extraigo las submatrices y especifico las cantidades conocidas
% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |    recuerde que siempre qc=0
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |    en este caso en particular fd=0

Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% desplazamientos para los gdls c = [1 2 3 10 11 12]
ac = [0; 0; 0; 0; 0; 0]; % desplazamientos conocidos

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);   % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd; % calculo fuerzas de equilibrio desconocidas

% armo los vectores de desplazamientos (a) y fuerzas (q)
a = zeros(ngdl,1);  q = zeros(ngdl,1);   % separo la memoria
a(c) = ac;   a(d) = ad; % desplazamientos 
q(c) = qd;  %q(d) = qc; % fuerzas nodales de equilibrio

fprintf('Desplazamientos de los nodos en coord. globales = \n'); a

%% imprimo las fuerzas internas en cada barra referidas a las coordenadas
%% globales
for e = 1:nb % para cada barra
   fprintf('\n\nFuerzas internas para barra %d en coord. globales = \n', e);
   qe_coord_glob = Ke{e}*a(idx(e,:)) - fe{e};
   disp(qe_coord_glob)
   
   fprintf('\nFuerzas internas para barra %d en coord. locales = \n', e);
   disp(T{e}*qe_coord_glob)
end;
