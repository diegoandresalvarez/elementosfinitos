% Ejemplo pag 282

clear
clc

EI = 1; % kN*m^2

%        b1         b2         b3      b4         b5
% v1 t2 ---- v3 t4 ---- v5 t6 ---- v7 ---- v8 t9 ---- v10 t11      >>>  GDLs
%   A          B          C         D                    E

ke  = cell(5,1);
idx = cell(5,1);

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

% Se ensambla la matriz de rigidez global   
ngdl = 11; nbar = 5;
K = zeros(ngdl);                    f = zeros(ngdl,1); 
                                    f(3) = -4; %kN
                                    f(4) = +4; %kN                                    
                                    f(8) = -2; %kN
                                    
for e = 1:nbar
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + ke{e};   
end

% Se separan los gdls conocidos de los desconocidos
c = [ 1 5 10 ]; d = setdiff(1:ngdl,c);

% Se resuelve el sistema  de ecuaciones
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = [ 0; 0; 0 ]; % desplazamientos conocidos en contorno

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc - Kdc*ac);      % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a  = zeros(ngdl,1);   a(c) = ac;  a(d) = ad; % desplazamientos
q  = zeros(ngdl,1);   q(c) = qd;             % fuerzas nodales equivalentes

a
q

for e = 1:nbar
   e, ke{e}*a(idx{e})
end