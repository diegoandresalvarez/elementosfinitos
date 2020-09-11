% Ejemplo pag 77

clear
clc

E = 200e9; % Pa
I = 50e-6; % m^4

% v1 t2               v3          v4 t5            v6 t7

% Barra 1
idx1 = [1 2 3];
EI = 2*E*I;
L = 4; %m
K1 = [ ...  
   3*EI/L^3,  3*EI/L^2,  -3*EI/L^3
   3*EI/L^2,    3*EI/L,  -3*EI/L^2
  -3*EI/L^3, -3*EI/L^2,   3*EI/L^3];
%{
clear, clc

syms syms V(x) M(x) t(x) v(x) EI L q

sol = dsolve(...       
      diff(V) == q,    diff(M) == V, ... % se definen las ecuaciones diferenciales
      diff(t) == M/EI, diff(v) == t, ... % 
      v(0) == 0, v(L) == 0,          ... % con sus respectivas condiciones
      t(0) == 0, M(L) == 0);             % de frontera.
      
% Se evaluan las cargas nodales equivalentes
disp('Y1 = '); disp(-subs(sol.V, x, 0));
disp('Y2 = '); disp(+subs(sol.V, x, L));

% Y los momentos nodales equivalentes
disp('M1 = '); disp(+subs(sol.M, x, 0));
disp('M2 = '); disp(-subs(sol.M, x, L));
%}
% Y1 = (5*L*q)/8; Y2 = (3*L*q)/8; M1 = (L^2*q)/8; M2 = 0;
q = -9e3; % N/m
f1 = [ (5*L*q)/8; (L^2*q)/8; (3*L*q)/8  ]; % Y1 M1 Y2

 % Barra 2
idx2 = [3 4 5]; 
EI = E*I;
L = 2; % m
K2 = [ ...
   3*EI/L^3, -3*EI/L^3,   3*EI/L^2
  -3*EI/L^3,  3*EI/L^3,  -3*EI/L^2
   3*EI/L^2, -3*EI/L^2,     3*EI/L];
f2 = [ 0; 0; 0];

% Barra 3
idx3 = [4 5 6 7];
EI = E*I;
L = 2; % m
K3 = [ ...
  12*EI/L^3,  6*EI/L^2, -12*EI/L^3,   6*EI/L^2   
   6*EI/L^2,    4*EI/L,  -6*EI/L^2,     2*EI/L
 -12*EI/L^3, -6*EI/L^2,  12*EI/L^3,  -6*EI/L^2
   6*EI/L^2,    2*EI/L,  -6*EI/L^2,     4*EI/L];



f3 = [ 0; 0; 0; 0];

% Se ensambla la matriz de rigidez global   
ngdl = 7;
K = zeros(ngdl);                    f = zeros(ngdl,1); 
                                    f(4) = -30e3; %N
K(idx1,idx1) = K1;                  f(idx1) = f1;
K(idx2,idx2) = K(idx2,idx2) + K2;   f(idx2) = f(idx2) + f2;
K(idx3,idx3) = K(idx3,idx3) + K3;   f(idx3) = f(idx3) + f3;

% Se separan los gdls conocidos de los desconocidos
c = [ 1 2 3 6 7]; d = [4 5];

% Se resuelve el sistema  de ecuaciones
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = [ 0; 0; 0; 0; 0 ]; % desplazamientos conocidos en contorno

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);        % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a  = zeros(ngdl,1);   a(c) = ac;  a(d) = ad; % desplazamientos
q  = zeros(ngdl,1);   q(c) = qd;             % fuerzas nodales equivalentes


a
q
