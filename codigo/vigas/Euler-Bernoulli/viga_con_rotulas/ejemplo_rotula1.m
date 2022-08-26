clear, clc

%% Se define una rigidez a la flexión cualquiera 
% (solo quiero calcular M y V)
EI = 1; % kN*m^2

%% 
% bx = barra
% vx = gdl desplazamiento vertical
% tx = gdl rotacion
%        b1         b2         b3      b4         b5
% v1 t2 ---- v3 t4 ---- v5 t6 ---- v7 ---- v8 t9 ---- v10 t11

%% Se separa memoria
ke  = cell(5,1);   % matrices de rigidez
idx = cell(5,1);   % gdl asociados a cada EF

%% Se especifican las matrices de rigidez de cada elemento de viga
% Barra 1
idx{1} = [1 2 3 4];
L  = 2; % m
ke{1} = [ ...
  12*EI/L^3,  6*EI/L^2, -12*EI/L^3,   6*EI/L^2   
   6*EI/L^2,    4*EI/L,  -6*EI/L^2,     2*EI/L
 -12*EI/L^3, -6*EI/L^2,  12*EI/L^3,  -6*EI/L^2
   6*EI/L^2,    2*EI/L,  -6*EI/L^2,     4*EI/L];

% Barra 2
idx{2} = [3 4 5 6];
L  = 2; % m
ke{2} = [ ...
  12*EI/L^3,  6*EI/L^2, -12*EI/L^3,   6*EI/L^2   
   6*EI/L^2,    4*EI/L,  -6*EI/L^2,     2*EI/L
 -12*EI/L^3, -6*EI/L^2,  12*EI/L^3,  -6*EI/L^2
   6*EI/L^2,    2*EI/L,  -6*EI/L^2,     4*EI/L];

% Barra 3
idx{3} = [5 6 7];
L  = 2; % m
ke{3} = [ ...  % rotula a la derecha
   3*EI/L^3,  3*EI/L^2,  -3*EI/L^3
   3*EI/L^2,    3*EI/L,  -3*EI/L^2
  -3*EI/L^3, -3*EI/L^2,   3*EI/L^3];

% Barra 4
idx{4} = [7 8 9];
L  = 1; % m
ke{4} = [ ...  % rotula a la izquierda
   3*EI/L^3, -3*EI/L^3,   3*EI/L^2
  -3*EI/L^3,  3*EI/L^3,  -3*EI/L^2
   3*EI/L^2, -3*EI/L^2,     3*EI/L];

% Barra 5
idx{5} = [8 9 10 11];
L  = 1; % m
ke{5} = [ ...
  12*EI/L^3,  6*EI/L^2, -12*EI/L^3,   6*EI/L^2   
   6*EI/L^2,    4*EI/L,  -6*EI/L^2,     2*EI/L
 -12*EI/L^3, -6*EI/L^2,  12*EI/L^3,  -6*EI/L^2
   6*EI/L^2,    2*EI/L,  -6*EI/L^2,     4*EI/L];

%% Se ensambla la matriz de rigidez global y el vector de fuerzas nodales equivalentes
ngdl = 11; nbar = 5;
K = zeros(ngdl);                    f = zeros(ngdl,1); 
                                    f(3) = -4; %N
                                    f(4) = +4; %N*m                                    
                                    f(8) = -2; %N                                                                        
                                    
for e = 1:nbar
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + ke{e};   
end

%% Se separan los gdls conocidos de los desconocidos
c = [ 1 5 10 ]; d = setdiff(1:ngdl,c);

%% Se resuelve el sistema  de ecuaciones
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = [ 0; 0; 0 ]; % desplazamientos conocidos en contorno

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);        % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a  = zeros(ngdl,1);   a(c) = ac;  a(d) = ad; % desplazamientos
q  = zeros(ngdl,1);   q(c) = qd;             % fuerzas nodales equivalentes

%% Se imprimen los resultados
a
q

format compact
for e = 1:nbar
   e
   qe = ke{e}*a(idx{e});
   switch e
       case {1, 2, 5}
           V = [  qe(1) -qe(3) ]
           M = [ -qe(2)  qe(4) ]
       case 3
           V = [  qe(1) -qe(3) ]
           M = [ -qe(2)  0     ]           
       case 4
           V = [  qe(1) -qe(2) ]
           M = [  0      qe(3) ]          
   end
end