% se definen algunas constantes que hacen el código más legible
NL1 = 1; NL2 = 2;

%% defino las variables
syms b E A L P          % define las variables simbolicas

long = [L;  L;  L/2]; % longitud de la barra

% LaG: local a global: matriz que relaciona nodos locales y globales
% (se lee la barra x va del nodo i al nodo j)
LaG = [1 3   % fila = barra
       2 3   % col1 = nodo global asociado a nodo local 1
       3 4]; % col2 = nodo global asociado a nodo local 2

k = E.*A./long;       % (k minuscula) rigidez de cada barra

% Fuerzas nodales equivalentes para cada barra
% fe = sym(zeros(2,3));
fe = [ b*L/2   b*L/2   0
       b*L/2   b*L/2   0 ];

%% Se prepara el vector de fuerzas nodales equivalentes global
f = sym(zeros(4,1));
f(3) = P/2;           % se agregan las cargas puntuales
f(4) = P;

%% ensamblo la matriz de rigidez global (K mayuscula)
K = sym(zeros(4));    % separa memoria para matriz de rigidez global K
for e = 1:3           % para cada una de las barras e = 1, 2 y 3
   idx = LaG(e,:);    % extrae indices de los nodos globales de la barra e
   Ke = k(e)*[1 -1; -1 1];            % matriz de rigidez local
   K(idx,idx) = K(idx,idx) + Ke;      % suma matriz de rigidez local
   f(idx)     = f(idx)     + fe(:,e); % suma vector de fuerzas nodales equivalentes local
end

%% Se imprime la matriz de rigidez
disp('Imprimamos la matriz de rigidez =')   
sympref('AbbreviateOutput', false);
pretty(K)
sympref('AbbreviateOutput','default');

%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
c = [1 2];    d = setdiff(1:4,c);

%% extraigo las submatrices y especifico las cantidades conocidas
% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |     recuerde siempre que qc=0
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |     en este caso en particular fd=0

Kcc = K(c,c); Kcd = K(c,d);
Kdc = K(d,c); Kdd = K(d,d);
fd  = f(c);   fc  = f(d);
ac = sym([0; 0]);

%% resuelvo el sistema de ecuaciones
% recuerde que \ es para resolver el sistema de ecuaciones eficientemente
ad = Kdd\(fc - Kdc*ac);  % = linsolve(Kdd, fc - Kdc*ac)
qd = Kcc*ac + Kcd*ad - fd;

%% formo los vectores de desplazamientos (a) y fuerzas (q)
a = sym(zeros(4,1)); q = sym(zeros(4,1));  % separo la memoria
a(c) = ac;           q(c) = qd;
a(d) = ad;         % q(d) = qc = 0;

%% calculo las cargas axiales en cada barra
fax = sym(zeros(3,1));
for e = 1:3
   fax(e) = k(e)*(a(LaG(e,NL2)) - a(LaG(e,NL1)));
end

%% imprimo los resultados
fprintf('\n\n a = \n'); pretty(a)
fprintf('\n\n q = \n'); pretty(q)
fprintf('\n\n fax = \n'); pretty(simplify(fax))
