% Ejemplo 11.3 Uribe Escamilla + soporte derecho inclinado 30 grados

%% Unidades en toneladas y cm

%% defino las variables
ang   = atan2d(300,400); % angulo especificado en grados

%barra    1       2       3       4       5
theta = [ ang     0       -ang    0       -90  ]; % angulo de inclinacion 
L     = [ 500     400     500     400      300 ]; % longitud barra
A     = [ 100     40      150     40       30  ]; % area barra

% LaG: local a global: matriz que relaciona nodos locales y globales
LaG = [1 3    % fila = barra
       1 4    % col1 = nodo global asociado a nodo local 1
       3 2    % col2 = nodo global asociado a nodo local 2
       4 2   
       3 4 ]; % (se lee la barra x va del nodo i al nodo j)

% gdl: grados de libertad
gdl = [ 1 2   % fila = nodo
        3 4   % col1 = gdl en direccion x
        5 6   % col2 = gdl en direccion y
        7 8 ]; 

E = 2040;     % ton/cm^2
k = E.*A./L;  % rigidez de cada barra

%% ensamblo la matriz de rigidez global
Kg = sparse(8,8); T = cell(5,1); % separo memoria
for e = 1:5  % para cada barra
   idx = [gdl(LaG(e,1),:) gdl(LaG(e,2),:)]; % saco los 4 gdls de la barra
   c = cosd(theta(e)); s = sind(theta(e));  % sin y cos de la inclinacion
   T{e} = [ c  s  0  0         % matriz de transformacion de coordenadas
           -s  c  0  0         % para la barra e
            0  0  c  s 
            0  0 -s  c ];
   Kloc = [ k(e)  0 -k(e)  0    % matriz de rigidez local expresada en el 
            0  0  0        0    % sistema de coordenadas locales para la                  
           -k(e)  0  k(e)  0    % barra e
            0     0  0     0 ];
   Kg(idx,idx) = Kg(idx,idx) + T{e}'*Kloc*T{e}; % sumo a K global
end

%% introduciendo los soportes inclinados
Tg = speye(8,8);
Tg(3:4,3:4) = [ cosd(30) sind(30); -sind(30) cosd(30) ];
K = Tg*Kg*Tg';

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
% desplazamientos para los gdls c = [1 2 4] en el sistema de coordenadas de
% los soportes inclinados
ac = [0; 0; 0]; 

% fuerzas en los gdls d = [3 5 6 7 8] en el sistema de coordenadas
% donde los grados de libertad son paralelos a los ejes x e y
fg  = zeros(8,1);
fcg = [0; 5*cosd(ang); 5*sind(ang); 0; -20]; %ton
fg(d) = fcg;

% convierto a sistema de coordenadas con soportes inclinados
f  = Tg*fg;
fc = f(d);

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc - Kdc*ac); % = linsolve(Kdd, fc - Kdc*ac)
qd = Kcc*ac + Kcd*ad;

% armo los vectores de desplazamientos (a) y fuerzas (q)
a = zeros(8,1); q = zeros(8,1);  % separo la memoria
a(c) = ac;       q(c) = qd;
a(d) = ad;     % q(d) = qc = 0

%% retorno las fuerzas y los desplazamientos en el sistema de coordenadas
%% donde los grados de libertad son paralelos a los ejes
qg = Tg'*q;
ag = Tg'*a;

%% calculo las fuerzas axiales (fax) en cada barra
fax = zeros(5,1);
for e = 1:5 % para cada barra
   idx = [gdl(LaG(e,1),:) gdl(LaG(e,2),:)]; % saco los 4 gdls de la barra e
   fax(e) = [-k(e) 0 k(e) 0]*T{e}*ag(idx);
end

%% imprimo los resultados
ag, qg, fax
