%% Ejemplo 11.3 Uribe Escamilla 

%% Unidades en toneladas y cm

%% se definen algunas constantes que hacen el código más legible
NL1 = 1; NL2 = 2;

%% defino las variables
ang   = atan2(300,400)*180/pi; % angulo especificado en grados

%barra    1       2       3       4       5
theta = [ ang     0       -ang    0       -90  ]; % angulo de inclinacion 
long  = [ 500     400     500     400      300 ]; % longitud barra
area  = [ 100     40      150     40       30  ]; % area barra

% LaG: local a global: matriz que relaciona nodos locales y globales
LaG = [1 3    % fila = barra
       1 4    % col1 = nodo global asociado a nodo local 1
       3 2    % col2 = nodo global asociado a nodo local 2
       4 2  
       3 4 ]; % (se lee la barra x va del nodo i al nodo j)

% gdl: grados de libertad
gdl = [1 2    % fila = nodo
       3 4    % col1 = gdl en direccion x
       5 6    % col2 = gdl en direccion y
       7 8]; 
      
% propiedades del material
E = 2040;           % ton/cm^2
k = E.*area./long;  % rigidez de cada barra

%% separo memoria
K = zeros(8); 
T = cell(5,1); 
idx = cell(5,1);

%% ensamblo la matriz de rigidez global
for e = 1:5  % para cada barra
   % saco los 4 gdls de la barra
   idx{e} = [gdl(LaG(e,NL1),:) gdl(LaG(e,NL2),:)]; 

   % matriz de transformacion de coordenadas para la barra e
   c = cosd(theta(e)); s = sind(theta(e));  % sin y cos de la inclinacion
   T{e} = [ c  s  0  0         
           -s  c  0  0
            0  0  c  s 
            0  0 -s  c ];

   % matriz de rigidez local expresada en el sistema de coordenadas locales para
   % la barra e
   Kloc = k(e)*[ 1  0 -1  0   
                 0  0  0  0
                -1  0  1  0
                 0  0  0  0 ];

   % sumo a K global
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + T{e}'*Kloc*T{e};
end

%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
c = [1 2 4];    d = [3 5 6 7 8];

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |    recuerde que siempre qc=0
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |    en este caso en particular fd=0

%% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d);
Kdc = K(d,c); Kdd = K(d,d);
% desplazamientos para los gdls c = [1 2 4]
ac = [0; 0; 0]; 

% fuerzas en los gdls d = [3 5 6 7 8]
fc = [0; 5*cosd(ang); 5*sind(ang); 0; -20]; %ton

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac); % = linsolve(Kdd, fc-Kdc*ac)
qd = Kcc*ac + Kcd*ad;

% armo los vectores de desplazamientos (a) y fuerzas (q)
a = zeros(8,1);  q = zeros(8,1);  % separo la memoria
a(c) = ac;       q(c) = qd;
a(d) = ad;     % q(d) = qc = 0

%% calculo las cargas axiales (N) en cada barra
N = zeros(5,1);
for e = 1:5 % para cada barra
   N(e) = k(e)*[-1 0 1 0]*T{e}*a(idx{e});
end

%% imprimo los resultados
a, q, N 
