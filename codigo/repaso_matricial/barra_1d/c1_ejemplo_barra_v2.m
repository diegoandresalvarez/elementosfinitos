% Ejemplo 1.1. Oñate

%% defino las variables
syms E A L P          % define las variables simbolicas

long = [L;  L;  L/2]; % longitud de la barra

% LaG: local a global: matriz que relaciona nodos locales y globales
% (se lee la barra x va del nodo i al nodo j)
LaG = [1 3   % fila = barra
       2 3   % col1 = nodo global asociado a nodo local 1
       3 4]; % col2 = nodo global asociado a nodo local 2

k = E.*A./long;       % (k minuscula) rigidez de cada barra
    
%% ensamblo la matriz de rigidez global (K mayuscula)
K = sym(zeros(4));    % separa memoria para matriz de rigidez global K
for e = 1:3           % para cada una de las barras e = 1, 2 y 3
   idx = LaG(e,:);    % extrae indices de los nodos globales de la barra e
   K(idx,idx) = K(idx,idx) + k(e)*[1 -1; -1 1]; % suma matriz de rigidez local
end;

%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
c = [1 2];    d = setdiff(1:4,c);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |     recuerde siempre que qc=0
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |     en este caso en particular fd=0

%% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); ac = sym([0; 0]);
Kdc = K(d,c); Kdd = K(d,d);                    fc = sym([0; P]);

%% resuelvo el sistema de ecuaciones
% recuerde que \ es para resolver el sistema de ecuaciones eficientemente
ad = Kdd\(fc-Kdc*ac);  % = linsolve(Kdd, fc-Kdc*ac)
qd = Kcc*ac + Kcd*ad;

%% formo los vectores de desplazamientos (a) y fuerzas (q)
a = sym(zeros(4,1)); q = sym(zeros(4,1));  % separo la memoria
a(c) = ac;           q(c) = qd;
a(d) = ad;         % q(d) = qc = 0;

%% calculo las cargas axiales en cada barra
N = sym(zeros(3,1));
for e = 1:3
   N(e) = k(e)*(a(LaG(e,2)) - a(LaG(e,1)));
end;

%% imprimo los resultados
fprintf('\n\n a = \n'); pretty(a)
fprintf('\n\n q = \n'); pretty(q)
fprintf('\n\n N = \n'); pretty(N)
